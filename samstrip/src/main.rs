use std::io::{self, BufRead, Write};

fn main() {
    if std::env::args().any(|arg| arg == "-h" || arg == "--help") {
        eprintln!("{}", HELP_MESSAGE);
        std::process::exit(1)
    }
    let mut stdout = io::stdout().lock();
    let mut stdin = io::stdin().lock();
    let mut line = String::new();
    let mut n_bytes = stdin.read_line(&mut line).unwrap();
    if n_bytes > 0 && !line.starts_with('@') {
        eprintln!(
            "Error: First SAM line did not start with a @, indicating a missing header. \
            \nDid you remember to pass in the whole SAM line? \
            If producing the SAM file using `samtools view`, remember the `-h` flag to include the header."
        );
        std::process::exit(1);
    }
    let mut buf = Vec::new();
    while n_bytes > 0 {
        if line.starts_with('@') {
            write!(&mut stdout, "{}", line).unwrap();
            line.clear();
            n_bytes = stdin.read_line(&mut line).unwrap();
            continue;
        }
        // This whole thing with `vec` vs `buf` is an elaborate workaround for
        // the borrowchecker. It needs to understand that, although `buf` is defined
        // outside this loop, and store references to the line in this loop,
        // it can't possibly hold any invalid references.
        // We do this with the awkward map unreachable at the end, which somehow
        // makes the compiler reset the lifetime of buf.
        let mut vec = buf;
        vec.extend(line.split('\t'));
        vec[9] = "*";
        vec[10] = "*";
        if let Some(i) = vec[11..].iter().position(|&s| s.starts_with("NM:i:")) {
            vec[11] = vec[i + 11];
            vec.truncate(12);
        } else {
            vec.truncate(11);
        }
        let mut whole_line = vec.join("\t");
        strip_line(&mut whole_line);
        writeln!(stdout, "{}", whole_line).unwrap();
        vec.clear();
        // See comment where vec is defined
        buf = vec.into_iter().map(|_| unreachable!()).collect();
        line.clear();
        n_bytes = stdin.read_line(&mut line).unwrap();
    }
}

fn strip_line(s: &mut String) {
    if s.ends_with('\n') {
        s.pop();
        if s.ends_with('\r') {
            s.pop();
        }
    }
}

const HELP_MESSAGE: &str = "samstrip

Reads a SAM file from stdin, and prints the equivalent stripped file to stdout.
A stripped file has the SEQ and QUAL fields removed, and all optional fields,
except the NM tag.
Barring any aligner-specific optional fields, a stripped SAM file contain the
same alignment information as a full file, but takes up less disk space.

Example usage: `samtools view -h file.bam | samstrip | samtools view -b - > stripped.bam`";
