use std::io::{self, BufRead, BufReader, BufWriter, Write};

struct FieldIterator<'a> {
    bytes: Option<&'a [u8]>,
}

impl<'a> FieldIterator<'a> {
    fn new(s: &'a [u8]) -> Self {
        Self { bytes: Some(s) }
    }
}

impl<'a> Iterator for FieldIterator<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(mem) = self.bytes {
            if let Some(next_pos) = memchr::memchr(b'\t', mem) {
                let (res, new) = mem.split_at(next_pos);
                // Safety: New is guaranteed to be nonempty by the semantics of split_at
                unsafe { self.bytes = Some(new.get_unchecked(1..new.len())) }
                Some(res)
            } else {
                self.bytes.take()
            }
        } else {
            None
        }
    }
}

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
    let mut stdout = BufWriter::new(io::stdout().lock());
    let mut stdin = BufReader::new(io::stdin().lock());
    match stdin.fill_buf().unwrap().first().map(|&b| b == b'@') {
        Some(true) => {}
        Some(false) => {
            eprintln!(
                "Error: First SAM line did not start with a @, indicating a missing header. \
                \nDid you remember to pass in the whole SAM file? \
                If producing the SAM file using `samtools view`, remember the `-h` flag to include the header."
            );
            std::process::exit(1);
        }
        None => std::process::exit(0),
    }

    let mut line: Vec<u8> = Vec::new();
    let mut n_bytes = stdin.read_until(b'\n', &mut line).unwrap();
    if line.first().map(|&b| b != b'@').unwrap_or(true) {
        eprintln!(
            "Error: First SAM line did not start with a @, indicating a missing header. \
            \nDid you remember to pass in the whole SAM file? \
            If producing the SAM file using `samtools view`, remember the `-h` flag to include the header."
        );
        std::process::exit(1);
    }
    while n_bytes > 0 {
        // TODO: We could speed this up by copying the header in a separate loop
        if line.first().map(|&b| b == b'@').unwrap_or(true) {
            stdout.write_all(&line).unwrap();
            line.clear();
            n_bytes = stdin.read_until(b'\n', &mut line).unwrap();
            continue;
        }
        // Write the line from the start to the 9th tab, as the first 9 fields
        // should be copied over unchanged
        let start = memchr::memchr_iter(b'\t', &line).nth(8).unwrap_or_else(|| {
            eprintln!("Error: In SAM alignment line, did not see all required fields");
            std::process::exit(1)
        });
        // The 10th and 11th field are SEQ and QUAL - we replace them with *
        stdout.write_all(&line[..=start]).unwrap();
        stdout.write_all(b"*\t*").unwrap();
        // Loop over the last fields - skip the first two SEQ and QUAL fields which we replaced above,
        // and then look for the NM:i field which should always be present. If we find it, copy it
        // to the out buffer.
        let fields = FieldIterator::new(&line[start + 1..]);
        let need_newline = if let Some(nm) = fields.skip(2).find(|&b| b.starts_with(b"NM:i")) {
            stdout.write_all(b"\t").unwrap();
            stdout.write_all(nm).unwrap();
            nm.last().map(|&b| b != b'\n').unwrap_or(true)
        } else {
            true
        };
        // The NM field which we copied might have ended in a newline, so check for existing newline
        // before adding one.
        if need_newline {
            stdout.write_all(b"\n").unwrap();
        }
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

Example usage:
Stripping an exiting BAM file:
`samtools view -h file.bam | samstrip | samtools view -b - > stripped.bam`

Stripping a BAM file while creating it:
`minimap2 -ax sr ref.fa fw.fq rv.fq | samstrip | samtools view -b - > file.bam`";
