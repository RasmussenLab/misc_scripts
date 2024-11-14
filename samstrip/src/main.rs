use std::io::{self, BufRead, Write};

fn main() {
    if std::env::args().any(|arg| arg == "-h" || arg == "--help") {
        eprintln!("{}", HELP_MESSAGE);
        std::process::exit(1)
    }
    let mut stdout = io::stdout().lock();
    let mut lines = io::stdin().lock().lines().map(Result::unwrap).peekable();
    if let Some(first_line) = lines.peek() {
        if !first_line.starts_with('@') {
            eprintln!(
                "Error: First SAM line did not start with a @, indicating a missing header. \
                \nDid you remember to pass in the whole SAM line? \
                If producing the SAM file using `samtools view`, remember the `-h` flag to include the header."
            );
            std::process::exit(1);
        }
    }
    for line in io::stdin().lock().lines().map(Result::unwrap) {
        if line.starts_with('@') {
            writeln!(&mut stdout, "{}", line).unwrap();
            continue
        }
        let mut buf: Vec<_> = line.split('\t').collect();
        buf[9] = "*";
        buf[10] = "*";
        if let Some(i) = buf[11..].iter().position(|&s| s.starts_with("NM:i:")) {
            buf[11] = buf[i + 11];
            buf.truncate(12);
        } else {
            buf.truncate(11);
        }
        writeln!(stdout, "{}", buf.join("\t")).unwrap();
    }
}

const HELP_MESSAGE: &'static str = "samstrip

Reads a SAM file from stdin, and prints the equivalent stripped file to stdout.
A stripped file has the SEQ and QUAL fields removed, and all optional fields,
except the NM tag.
Barring any aligner-specific optional fields, a stripped SAM file contain the
same alignment information as a full file, but takes up less disk space.

Example usage: `samtools view -h file.bam | samstrip | samtools view -b - > stripped.bam`";