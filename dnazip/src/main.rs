use anyhow::{self, Context};
use clap::Parser;

use crossbeam_channel::{self, Receiver, RecvError, TryRecvError};
use flate2::{bufread::GzEncoder, Compression};
use size;
use std::io::{stderr, ErrorKind, Write};
use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::{Path, PathBuf},
    thread,
};
use walkdir::WalkDir;

/// Gzip compresses all FAST{Q,A} files found recursively in the given directory.
/// Does not follow symlinks.
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Directory to start from
    start: PathBuf,

    /// Print the paths that would be compressed; do not compress
    #[arg(short, long)]
    dry_run: bool,

    /// Print paths that are being compressed
    #[arg(short, long)]
    verbose: bool,

    /// Number of additional threads to use for compression [0]
    #[arg(short, long, default_value_t = 0)]
    threads: u8,
}

const FASTA_EXTENSIONS: [&str; 4] = ["fna", "fasta", "fa", "faa"];
const FASTQ_EXTENSIONS: [&str; 2] = ["fq", "fastq"];

fn is_fasta(p: &Path) -> bool {
    p.extension()
        .is_some_and(|e| e.to_str().is_some_and(|s| FASTA_EXTENSIONS.contains(&s)))
}

fn is_fastq(p: &Path) -> bool {
    p.extension()
        .is_some_and(|e| e.to_str().is_some_and(|s| FASTQ_EXTENSIONS.contains(&s)))
}

fn write_path(path: &Path, prefix: Option<&str>) {
    let mut v: Vec<u8> = Vec::new();
    if let Some(s) = prefix {
        v.write_all(s.as_bytes()).unwrap();
    }
    v.write_all(path.as_os_str().as_encoded_bytes()).unwrap();
    v.write_all(&[b'\n']).unwrap();
    stderr().write_all(&v).unwrap()

}

fn compress(path: &Path, dry_run: bool, verbose: bool) -> anyhow::Result<()> {
    if dry_run {
        write_path(path, Some("Would compress: "));
        return Ok(())
    }
    let mut p = path.as_os_str().to_owned();
    p.push(".gz");
    let mut dst = BufWriter::new(
        File::create(&p).with_context(|| format!("Could not create gzipped file {:?}", p))?,
    );
    let mut new = GzEncoder::new(
        BufReader::new(
            File::open(path).with_context(|| format!("Could not open file: {:?}", path))?,
        ),
        Compression::default(),
    );
    std::io::copy(&mut new, &mut dst).context("Error when copying file to gzip wrier")?;
    std::fs::remove_file(path).with_context(|| format!("Could not remove file {:?}", path))?;
    if verbose {
        write_path(path, Some("Compressed: "))
    }
    Ok(())
}

fn read_channel(reciever: Receiver<PathBuf>, dry_run: bool, verbose: bool) {
    loop {
        match reciever.recv() {
            Err(RecvError) => return,
            Ok(path) => compress(&path, dry_run, verbose).unwrap(),
        }
    }
}

fn main() {
    let args = Cli::parse();
    let mut n_files = 0;
    let mut n_bytes = 0;
    // We cap the channel to hold a number of entries which is low enough that it won't
    let (sender, reciever) = crossbeam_channel::unbounded::<PathBuf>();
    let handles: Vec<_> = (0..args.threads)
        .map(|_| {
            let rec = reciever.clone();
            thread::spawn(move || read_channel(rec, args.dry_run, args.verbose))
        })
        .collect();
    for maybe_entry in WalkDir::new(args.start) {
        let handled_entry = match maybe_entry {
            Ok(e) => Some(Ok(e)),
            Err(err) => {
                let path = err.path().unwrap_or(Path::new("")).display();
                if let Some(inner) = err.io_error() {
                    match inner.kind() {
                        ErrorKind::PermissionDenied => {
                            eprintln!("Permission denied: {}", path);
                            None
                        }
                        _ => Some(Err(err)),
                    }
                } else {
                    Some(Err(err))
                }
            }
        };
        let entry = if let Some(res) = handled_entry {
            res.unwrap()
        } else {
            continue;
        };
        if entry.file_type().is_file() && !entry.path_is_symlink() {
            let path = entry.path();
            if is_fasta(path) || is_fastq(path) {
                sender.send(path.to_owned()).unwrap();
                n_files += 1;
                n_bytes += entry.metadata().unwrap().len()
            }
        } else {
            continue;
        }
        // If there are no dedicated readers, we use the main thread to compress an entry.
        // This way the main thread never has to wait for the worker threads.
        if args.threads == 0 {
            match reciever.try_recv() {
                Ok(p) => compress(&p, args.dry_run, args.verbose).unwrap(),
                Err(TryRecvError::Disconnected) => unreachable!(),
                // Below can also never happen, but no big deal if it does
                Err(TryRecvError::Empty) => (),
            }
        }
    }

    // This signals to the worker threads that they should exit,
    // once they run out of paths to process
    drop(sender);

    // Turn the main thread into a worker thread to help with the last paths.
    read_channel(reciever, args.dry_run, args.verbose);

    // Make sure all the workers exited
    for handle in handles {
        handle.join().unwrap()
    }
    if args.dry_run {
        eprintln!(
            "Would compress {} files, {}",
            n_files,
            size::Size::from_bytes(n_bytes)
        );
    } else {
        eprintln!(
            "Compressed {} files, {}",
            n_files,
            size::Size::from_bytes(n_bytes)
        );
    }
}
