//! Rust port of biofast's `fqcnt_c1_kseq.c` and the FASTA/FASTQ parser of
//! `kseq.h`: count the records, bases, and quality values of a file.
//!
//! `BufReader` plays the part of `kstream_t`; the parser is written against
//! `BufRead` and follows `kseq_read`, including multi-line records.

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, ErrorKind};
use std::process;

const BUFFER_SIZE: usize = 16384;

/// `kseq_t`: the most recently read record and the stream it came from.
struct Kseq<R> {
    name: Vec<u8>,
    comment: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    /// The `>` or `@` of the next record, when reading a sequence consumed it.
    last_char: Option<u8>,
    reader: R,
}

/// `isspace` in the C locale.
fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t'..=b'\r')
}

fn truncated_quality() -> io::Error {
    io::Error::new(ErrorKind::InvalidData, "truncated quality string")
}

impl<R: BufRead> Kseq<R> {
    fn new(reader: R) -> Self {
        Kseq {
            name: Vec::new(),
            comment: Vec::new(),
            seq: Vec::with_capacity(256),
            qual: Vec::new(),
            last_char: None,
            reader,
        }
    }

    /// `ks_getc`.
    fn getc(&mut self) -> io::Result<Option<u8>> {
        let byte = self.reader.fill_buf()?.first().copied();
        if byte.is_some() {
            self.reader.consume(1);
        }
        Ok(byte)
    }

    /// `ks_getuntil2` with `KS_SEP_SPACE`: append the bytes before the next
    /// whitespace to `text` and return that whitespace, or `None` at the end
    /// of the input.
    fn read_word(reader: &mut R, text: &mut Vec<u8>) -> io::Result<Option<u8>> {
        loop {
            let buffer = reader.fill_buf()?;
            if buffer.is_empty() {
                return Ok(None);
            }
            match buffer.iter().position(|&byte| is_space(byte)) {
                Some(i) => {
                    let delimiter = buffer[i];
                    text.extend_from_slice(&buffer[..i]);
                    reader.consume(i + 1);
                    return Ok(Some(delimiter));
                }
                None => {
                    let length = buffer.len();
                    text.extend_from_slice(buffer);
                    reader.consume(length);
                }
            }
        }
    }

    /// `ks_getuntil2` with `KS_SEP_LINE`: append the rest of the line to
    /// `text`, without its terminator. Returns whether the input had anything
    /// left to read.
    fn read_line(reader: &mut R, text: &mut Vec<u8>) -> io::Result<bool> {
        if reader.read_until(b'\n', text)? == 0 {
            return Ok(false);
        }
        if text.last() == Some(&b'\n') {
            text.pop();
        }
        if text.len() > 1 && text.last() == Some(&b'\r') {
            text.pop();
        }
        Ok(true)
    }

    /// `kseq_read`: read the next record into `self`. Returns `false` at the
    /// end of the input and an `InvalidData` error for a FASTQ record whose
    /// quality string is missing or not as long as its sequence.
    fn read(&mut self) -> io::Result<bool> {
        if self.last_char.is_none() {
            // Jump to the next header line.
            loop {
                match self.getc()? {
                    None => return Ok(false),
                    Some(byte @ (b'>' | b'@')) => {
                        self.last_char = Some(byte);
                        break;
                    }
                    Some(_) => {}
                }
            }
        }
        self.name.clear();
        self.comment.clear();
        self.seq.clear();
        self.qual.clear();

        let delimiter = Self::read_word(&mut self.reader, &mut self.name)?;
        if delimiter.is_none() && self.name.is_empty() {
            return Ok(false);
        }
        if delimiter != Some(b'\n') {
            Self::read_line(&mut self.reader, &mut self.comment)?;
        }

        let next = loop {
            match self.getc()? {
                Some(b'\n') => {} // Skip empty lines.
                Some(byte) if !matches!(byte, b'>' | b'+' | b'@') => {
                    self.seq.push(byte);
                    Self::read_line(&mut self.reader, &mut self.seq)?;
                }
                next => break next,
            }
        };
        if let Some(b'>' | b'@') = next {
            self.last_char = next; // The next record's header has begun.
        }
        if next != Some(b'+') {
            return Ok(true); // FASTA
        }

        // Skip the rest of the '+' line.
        loop {
            match self.getc()? {
                None => return Err(truncated_quality()),
                Some(b'\n') => break,
                Some(_) => {}
            }
        }
        while Self::read_line(&mut self.reader, &mut self.qual)? && self.qual.len() < self.seq.len() {}
        self.last_char = None; // The next header line is still ahead.
        if self.seq.len() != self.qual.len() {
            return Err(truncated_quality());
        }
        Ok(true)
    }
}

fn main() {
    let Some(path) = env::args().nth(1) else {
        print!("Usage: fqcnt <in.fq>");
        return;
    };
    let file = File::open(&path).unwrap_or_else(|error| {
        eprintln!("fqcnt: {path}: {error}");
        process::exit(1);
    });
    let mut record = Kseq::new(BufReader::with_capacity(BUFFER_SIZE, file));

    let (mut count, mut seq_length, mut qual_length) = (0u64, 0usize, 0usize);
    let result = loop {
        match record.read() {
            Ok(true) => {
                count += 1;
                seq_length += record.seq.len();
                qual_length += record.qual.len();
            }
            Ok(false) => break Ok(()),
            Err(error) => break Err(error),
        }
    };
    println!("{count}\t{seq_length}\t{qual_length}");
    if result.is_err() {
        eprintln!("ERROR: malformated FASTX");
    }
}
