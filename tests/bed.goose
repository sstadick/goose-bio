import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";
import bio from "..";
import bio.bed;

fn main() {
    let file = "./tests/example.bed";
    let chr_map = bed::ChrLookup {};
    let ret = bed::read_bedfile(file, chr_map) { record, rest =>
        print(record);
        print("Rest: ", rest);
    };
    assert(ret);
}
