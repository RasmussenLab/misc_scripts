const HELP: &str = "bbformat: \
Convert a Vamb .tsv output binning file to CAMI Bioboxes binning format, \
which is used as input to AMBER. \
Reads from stdin, prints to stdout.
Example use: cat vamb_clusters.tsv | bbformat > amber.bb";

// To compile: Install Rust 1.73 or newer.
// Compile with `rustc -C strip=debuginfo -C opt-level=s bbformat.rs`

use std::io::Write;

fn print_help_and_exit() -> ! {
    eprintln!("{}", HELP);
    std::process::exit(1)
}

fn print_usage_and_exit() -> ! {
    eprintln!("Usage: cat input.tsv | bbformat > out.bb.\nRun bbformat -h for help");
    std::process::exit(1)
}

fn main() {
    let mut too_many_args = false;
    for (argno, arg) in std::env::args().enumerate() {
        if arg == "--help" || arg == "-h" {
            print_help_and_exit()
        }
        too_many_args |= argno > 0;
    }
    if too_many_args {
        print_usage_and_exit()
    }
    // Note that the output format MUST be UTF8 per the specs, and so must the input also be.
    let s = std::io::read_to_string(std::io::stdin()).expect("Could not read input file as UTF-8");
    let mut lines = s.trim().lines().peekable();
    // Skip the header if it exists
    lines.next_if_eq(&"clustername\tcontigname");
    let mut stdout = std::io::stdout().lock();
    stdout
        .write_all(b"@Version:0.9.1\n@SampleID:all\n\n@@SEQUENCEID\tBINID\n")
        .expect("Unable to write header");
    for line in lines {
        let (cluster, contig) = line
            .split_once('\t')
            .expect("Expected a tab character on line");
        if contig.as_bytes().contains(&b'\t') {
            panic!("Input line has more than two tab-separated fields")
        }
        writeln!(stdout, "{}\t{}", contig, cluster).expect("Unable to write to output file");
    }
}
