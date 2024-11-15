use std::io::{self, BufRead, Write};

fn main() {
    if let Some(arg) = std::env::args().nth(1) {
        if arg == "-h" || arg == "--help" {
            eprintln!("{}", HELP_MESSAGE);
            std::process::exit(0)
        } else if arg == "-V" || arg == "--version" {
            eprintln!("samstrip version {}", VERSION);
            std::process::exit(0)
        } else {
            eprintln!("Unrecognized argument: {}.\nRun samstrip -h for help.", arg);
            std::process::exit(1)
        }
    }
    let mut stdout = io::stdout().lock();
    let mut stdin = io::stdin().lock();
    let mut line: Vec<u8> = Vec::new();
    let mut to_write: Vec<u8> = Vec::new();
    let mut n_bytes = stdin.read_until(b'\n', &mut line).unwrap();
    if line.first().map(|&b| b != b'@').unwrap_or(true) {
        eprintln!(
            "Error: First SAM line did not start with a @, indicating a missing header. \
            \nDid you remember to pass in the whole SAM line? \
            If producing the SAM file using `samtools view`, remember the `-h` flag to include the header."
        );
        std::process::exit(1);
    }
    while n_bytes > 0 {
        if line.first().map(|&b| b == b'@').unwrap_or(true) {
            stdout.write_all(&line).unwrap();
            line.clear();
            n_bytes = stdin.read_until(b'\n', &mut line).unwrap();
            continue;
        }
        let mut fields = line.split(|&b| b == b'\t');
        if let Some(field) = fields.next() {
            to_write.extend_from_slice(field);
        }
        for field in fields.by_ref().take(8) {
            to_write.push(b'\t');
            to_write.extend_from_slice(field);
        }
        to_write.extend_from_slice(b"\t*\t*\t");
        if let Some(nm) = fields.skip(2).filter(|s| s.starts_with(b"NM")).next() {
            to_write.extend_from_slice(nm);
        }
        if !to_write.last().map(|&b| b == b'\n').unwrap_or(false) {
            to_write.push(b'\n')
        }
        stdout.write_all(&to_write).unwrap();
        to_write.clear();
        line.clear();
        n_bytes = stdin.read_until(b'\n', &mut line).unwrap();
    }
}

const VERSION: &str = env!("CARGO_PKG_VERSION");

const HELP_MESSAGE: &str = "samstrip

Reads a SAM file from stdin, and prints the equivalent stripped file to stdout.
A stripped file has the SEQ and QUAL fields removed, and all optional fields,
except the NM tag.
Barring any aligner-specific optional fields, a stripped SAM file contain the
same alignment information as a full file, but takes up less disk space.

Options:
    -h --help: Print help
    -V --version: Print version

Example usage: `samtools view -h file.bam | samstrip | samtools view -b - > stripped.bam`";
