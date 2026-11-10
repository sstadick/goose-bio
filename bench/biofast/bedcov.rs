//! Rust port of biofast's `bedcov_c1_cgr.c` and the parts of `cgranges.c` it
//! uses: load one BED file into an implicit interval tree, then print, for
//! every interval of a second BED file, how many loaded intervals overlap it
//! and how many of its bases they cover.
//!
//! The port is safe Rust on the standard library. It follows the C line by
//! line where that reads naturally and uses `Vec`, `HashMap`, `BufReader`,
//! and `BufWriter` where C manages memory and buffers by hand.

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::process;

const BUFFER_SIZE: usize = 0x10000;

// ---------------------------------------------------------------------------
// Radix sort
// ---------------------------------------------------------------------------

const MIN_SIZE: usize = 64;
const RADIX_BITS: u32 = 8;
const BUCKET_COUNT: usize = 1 << RADIX_BITS;

#[inline(always)]
fn digit(interval: Interval, shift: u32) -> usize {
    ((interval.x >> shift) & (BUCKET_COUNT as u64 - 1)) as usize
}

fn insertion_sort(slice: &mut [Interval]) {
    for i in 1..slice.len() {
        if slice[i].x < slice[i - 1].x {
            let temporary = slice[i];
            let mut j = i;
            while j > 0 && temporary.x < slice[j - 1].x {
                slice[j] = slice[j - 1];
                j -= 1;
            }
            slice[j] = temporary;
        }
    }
}

/// `rs_sort_cr_intv`, with bucket cursors and ends as two tables of indices
/// (see `bench/RADIX_SORT.md` for why this layout matches C's pointer pairs).
fn radix_sort_inner(slice: &mut [Interval], shift: u32) {
    let mut begins = [0usize; BUCKET_COUNT];
    let mut ends = [0usize; BUCKET_COUNT];
    for interval in slice.iter() {
        ends[digit(*interval, shift)] += 1;
    }
    let mut previous_end = 0;
    for bucket in 0..BUCKET_COUNT {
        begins[bucket] = previous_end;
        previous_end += ends[bucket];
        ends[bucket] = previous_end;
    }

    for bucket in 0..BUCKET_COUNT {
        while begins[bucket] != ends[bucket] {
            let origin = begins[bucket];
            let mut temporary = slice[origin];
            let mut target = digit(temporary, shift);
            while target != bucket {
                let destination = &mut begins[target];
                temporary = std::mem::replace(&mut slice[*destination], temporary);
                *destination += 1;
                target = digit(temporary, shift);
            }
            slice[origin] = temporary;
            begins[bucket] = origin + 1;
        }
    }

    if shift != 0 {
        let next_shift = shift.saturating_sub(RADIX_BITS);
        let mut begin = 0;
        for end in ends {
            let part = &mut slice[begin..end];
            if part.len() > MIN_SIZE {
                radix_sort_inner(part, next_shift);
            } else if part.len() > 1 {
                insertion_sort(part);
            }
            begin = end;
        }
    }
}

fn radix_sort(slice: &mut [Interval]) {
    if slice.len() <= MIN_SIZE {
        insertion_sort(slice);
    } else {
        radix_sort_inner(slice, u64::BITS - RADIX_BITS);
    }
}

// ---------------------------------------------------------------------------
// cgranges
// ---------------------------------------------------------------------------

/// `cr_intv_t`. Before indexing, `x` is `contig << 32 | start`; afterwards it
/// is `start << 32 | end` and `y` is the largest end in the interval's subtree.
#[derive(Clone, Copy)]
struct Interval {
    x: u64,
    y: i32,
    #[allow(dead_code)] // Carried, as in C, so a record is 16 bytes.
    label: i32,
}

impl Interval {
    #[inline(always)]
    fn start(&self) -> i32 {
        (self.x >> 32) as i32
    }

    #[inline(always)]
    fn end(&self) -> i32 {
        self.x as i32
    }
}

/// `cr_ctg_t` without the name, which the map owns, and the unused length.
#[derive(Clone, Copy)]
struct Contig {
    /// Height of the contig's tree, or -1 for a contig without intervals.
    root_k: i32,
    /// The contig's intervals are `intervals[offset..offset + count]`.
    offset: usize,
    count: usize,
}

#[derive(Default)]
struct Cgranges {
    intervals: Vec<Interval>,
    contigs: Vec<Contig>,
    contig_ids: HashMap<Vec<u8>, usize>,
}

/// One node of the traversal in `overlap`: the subtree rooted at index `x`
/// with height `k`, and whether its left child has been handled.
#[derive(Clone, Copy, Default)]
struct StackEntry {
    x: usize,
    k: i32,
    left_done: bool,
}

impl Cgranges {
    /// `cr_add_ctg`.
    fn add_contig(&mut self, name: &[u8]) -> usize {
        if let Some(&id) = self.contig_ids.get(name) {
            return id;
        }
        let id = self.contigs.len();
        self.contigs.push(Contig { root_k: -1, offset: 0, count: 0 });
        self.contig_ids.insert(name.to_vec(), id);
        id
    }

    /// `cr_add`.
    fn add(&mut self, contig: &[u8], start: i32, end: i32, label: i32) {
        if start > end {
            return;
        }
        let id = self.add_contig(contig);
        self.intervals.push(Interval { x: (id as u64) << 32 | start as u32 as u64, y: end, label });
    }

    /// `cr_index`: sort by contig and start, then build each contig's tree.
    fn index(&mut self) {
        // cr_index_prepare
        let intervals = &mut self.intervals;
        if !intervals.is_sorted_by_key(|interval| interval.x) {
            radix_sort(intervals);
        }
        let mut first = 0;
        for i in 1..=intervals.len() {
            if i == intervals.len() || intervals[i].x >> 32 != intervals[first].x >> 32 {
                let contig = &mut self.contigs[(intervals[first].x >> 32) as usize];
                contig.offset = first;
                contig.count = i - first;
                first = i;
            }
        }
        for interval in intervals.iter_mut() {
            interval.x = interval.x << 32 | interval.y as u32 as u64;
            interval.y = 0;
        }

        for contig in &mut self.contigs {
            let range = contig.offset..contig.offset + contig.count;
            contig.root_k = index_contig(&mut self.intervals[range]);
        }
    }

    /// `cr_overlap`: replace `hits` with the indices in `intervals`, in start
    /// order, of the intervals of `contig` that overlap `start..end`.
    fn overlap(&self, contig: &[u8], start: i32, end: i32, hits: &mut Vec<usize>) {
        hits.clear();
        let Some(&id) = self.contig_ids.get(contig) else { return };
        let contig = self.contigs[id];
        let tree = &self.intervals[contig.offset..contig.offset + contig.count];
        if contig.root_k < 0 {
            return;
        }

        let mut stack = [StackEntry::default(); 64];
        let mut depth = 1;
        stack[0] = StackEntry { x: (1 << contig.root_k) - 1, k: contig.root_k, left_done: false };
        while depth > 0 {
            depth -= 1;
            let node = stack[depth];
            if node.k <= 3 {
                // The subtree has at most (1 << (k + 1)) - 1 intervals: scan it.
                let first = node.x >> node.k << node.k;
                let last = (first + (1 << (node.k + 1)) - 1).min(tree.len());
                for i in first..last {
                    if tree[i].start() >= end {
                        break;
                    }
                    if start < tree[i].end() {
                        hits.push(contig.offset + i);
                    }
                }
            } else if !node.left_done {
                let left = node.x - (1 << (node.k - 1));
                stack[depth] = StackEntry { left_done: true, ..node };
                depth += 1;
                if left >= tree.len() || tree[left].y > start {
                    stack[depth] = StackEntry { x: left, k: node.k - 1, left_done: false };
                    depth += 1;
                }
            } else if node.x < tree.len() && tree[node.x].start() < end {
                if start < tree[node.x].end() {
                    hits.push(contig.offset + node.x);
                }
                let right = node.x + (1 << (node.k - 1));
                stack[depth] = StackEntry { x: right, k: node.k - 1, left_done: false };
                depth += 1;
            }
        }
    }
}

/// `cr_index1`: store in every `y` the largest end in that interval's subtree
/// and return the tree's height.
fn index_contig(tree: &mut [Interval]) -> i32 {
    if tree.is_empty() {
        return -1;
    }
    let count = tree.len();
    let mut last_i = 0;
    let mut last = 0;
    for i in (0..count).step_by(2) {
        last_i = i;
        last = tree[i].end();
        tree[i].y = last;
    }
    let mut k = 1;
    while 1 << k <= count {
        let x = 1 << (k - 1);
        for i in ((x << 1) - 1..count).step_by(x << 2) {
            let left = tree[i - x].y;
            let right = if i + x < count { tree[i + x].y } else { last };
            tree[i].y = tree[i].end().max(left).max(right);
        }
        last_i = if last_i >> k & 1 != 0 { last_i - x } else { last_i + x };
        if last_i < count && tree[last_i].y > last {
            last = tree[last_i].y;
        }
        k += 1;
    }
    k - 1
}

// ---------------------------------------------------------------------------
// bedcov
// ---------------------------------------------------------------------------

/// C's `atol` narrowed to 32 bits: optional blanks and sign, then digits.
fn atoi(text: &[u8]) -> i32 {
    let mut rest = text;
    while let [b' ' | b'\t'..=b'\r', tail @ ..] = rest {
        rest = tail;
    }
    let negative = match rest {
        [sign @ (b'+' | b'-'), tail @ ..] => {
            rest = tail;
            *sign == b'-'
        }
        _ => false,
    };
    let mut value = 0i64;
    for &byte in rest.iter().take_while(|byte| byte.is_ascii_digit()) {
        value = value.wrapping_mul(10).wrapping_add((byte - b'0') as i64);
    }
    (if negative { value.wrapping_neg() } else { value }) as i32
}

/// `parse_bed`: the contig, start, and end of a line with at least three
/// tab-separated fields.
fn parse_bed(line: &[u8]) -> Option<(&[u8], i32, i32)> {
    let mut fields = line.split(|&byte| byte == b'\t');
    let contig = fields.next()?;
    let start = atoi(fields.next()?);
    let end = atoi(fields.next()?);
    Some((contig, start, end))
}

/// Call `visit` with each line of the file, without its line terminator.
fn for_each_line(path: &str, mut visit: impl FnMut(&[u8]) -> io::Result<()>) -> io::Result<()> {
    let mut reader = BufReader::with_capacity(BUFFER_SIZE, File::open(path)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            return Ok(());
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.len() > 1 && line.last() == Some(&b'\r') {
            line.pop();
        }
        visit(&line)?;
    }
}

fn read_bed(path: &str) -> io::Result<Cgranges> {
    let mut ranges = Cgranges::default();
    let mut label = 0;
    for_each_line(path, |line| {
        if let Some((contig, start, end)) = parse_bed(line) {
            ranges.add(contig, start, end, label);
            label += 1;
        }
        Ok(())
    })?;
    Ok(ranges)
}

fn run(loaded: &str, streamed: &str) -> io::Result<()> {
    let mut ranges = read_bed(loaded)?;
    ranges.index();

    let mut out = BufWriter::with_capacity(BUFFER_SIZE, io::stdout().lock());
    let mut hits = Vec::new();
    for_each_line(streamed, |line| {
        let Some((contig, start1, end1)) = parse_bed(line) else { return Ok(()) };
        ranges.overlap(contig, start1, end1, &mut hits);
        let (mut coverage, mut covered_start, mut covered_end) = (0i64, 0i64, 0i64);
        for &hit in &hits {
            let interval = ranges.intervals[hit];
            let start0 = interval.start().max(start1) as i64;
            let end0 = interval.end().min(end1) as i64;
            if start0 > covered_end {
                coverage += covered_end - covered_start;
                covered_start = start0;
                covered_end = end0;
            } else {
                covered_end = covered_end.max(end0);
            }
        }
        coverage += covered_end - covered_start;
        out.write_all(contig)?;
        writeln!(out, "\t{start1}\t{end1}\t{}\t{coverage}", hits.len())
    })?;
    out.flush()
}

fn main() {
    let arguments: Vec<String> = env::args().collect();
    let [_, loaded, streamed] = arguments.as_slice() else {
        println!("Usage: bedcov <loaded.bed> <streamed.bed>");
        return;
    };
    if let Err(error) = run(loaded, streamed) {
        eprintln!("bedcov: {error}");
        process::exit(1);
    }
}
