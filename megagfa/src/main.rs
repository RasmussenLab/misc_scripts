use anyhow::{self, bail, Context, Result};
use bio::io::fasta::Reader;
use clap::Parser;
use regex::Regex;
use smallvec::SmallVec;
use std::{
    collections::HashMap,
    io::{stdin, stdout, BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

fn exitwith(s: &str) {
    eprintln!("{}", s);
    std::process::exit(1)
}

fn main() -> Result<()> {
    let args = Cli::parse();
    if args.k == 0 {
        exitwith("The value of K cannot be zero")
    }
    let input: Box<dyn BufRead> = if let Some(p) = args.i {
        if !p.is_file() {
            exitwith(&format!(
                "Not an existing file: \"{}\"",
                p.to_string_lossy()
            ));
        }
        p.file_name().and_then(|f| f.to_str()).and_then(|s| {
            let re = Regex::new(r"^k(\d+)\.contigs\.fa$").unwrap();
            re.captures(s).map(|c| (c, s))
        }).map(|(caps, s)| {
            let u = caps[1].parse::<u8>().ok()?;
            if u != args.k {
                exitwith(&format!("ERROR: K value passed with -k is {}, but given file is {} with different K value.", args.k, s))
            }
            None::<()>
        });
        Box::new(BufReader::new(
            std::fs::File::open(p.clone()).with_context(|| {
                format!("Could not open input file \"{}\"", p.to_string_lossy())
            })?,
        ))
    } else {
        Box::new(BufReader::new(stdin().lock()))
    };
    let (identifiers, edges) = build_graph(input, args.k, args.min_contig_length as usize)?;
    write_gfa(&identifiers, &edges)?;
    Ok(())
}

#[derive(Clone, Copy)]
struct KmerData(u32);

impl KmerData {
    fn try_new(index: usize, rc: bool) -> Result<Self> {
        let x: u32 = index
            .try_into()
            .ok()
            .and_then(|u| if u > 0x7fffffff { None } else { Some(u) })
            .context("Can only hande 2147483647 FASTA records in one file")?;
        Ok(Self(x + 0x80000000 * (rc as u32)))
    }

    fn is_rc(&self) -> bool {
        self.0 > 0x7fffffff
    }

    fn index(&self) -> usize {
        (self.0 & 0x7fffffff) as usize
    }
}

struct Edge {
    from: KmerData,
    to: KmerData,
}

const fn make_lut() -> [u8; 256] {
    let mut lut = [0xff; 256];
    let mut i: u8 = 0;
    while i < 127 {
        i += 1;
        lut[i as usize] = match i {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            _ => 0xff,
        }
    }
    lut
}

const LUT: [u8; 256] = make_lut();

// TODO: Custom 2-bit packed slice type instead
#[derive(PartialEq, Eq, Hash)]
struct Kmer(Vec<u8>);

impl Kmer {
    fn make_pair(bytes: &[u8], into: &mut [u8]) -> Option<(Self, Self)> {
        translate(bytes, into)?;
        let fw = compact(into);
        reverse_complement(into);
        let rv = compact(into);
        Some((Kmer(fw), Kmer(rv)))
    }
}

fn translate(s: &[u8], into: &mut [u8]) -> Option<()> {
    for (byte, to) in s.iter().zip(into.iter_mut()) {
        unsafe {
            let b = *LUT.get_unchecked(*byte as usize);
            if b == 0xff {
                return None
            } else {
                *to = b
            }
        }
    }
    Some(())
}

fn reverse_complement(x: &mut [u8]) {
    x.reverse();
    for i in x.iter_mut() {
        *i ^= 3
    }
}

fn compact(x: &[u8]) -> Vec<u8> {
    let chunks = x.chunks_exact(4);
    let r = chunks.remainder();
    let len = chunks.len() + !r.is_empty() as usize;
    let mut v: Vec<u8> = vec![0; len];
    for (chunk, i) in chunks.zip(v.iter_mut()) {
        *i = compact_element(chunk)
    }
    if !r.is_empty() {
        *v.last_mut().unwrap() = compact_element(r);
    };
    v
}

fn compact_element(s: &[u8]) -> u8 {
    s.iter().fold(0, |old, new| (old << 2) | new)
}

fn build_graph(
    input: impl BufRead,
    k: u8,
    min_contig_length: usize,
) -> Result<(Vec<Option<String>>, Vec<Edge>)> {
    let id_regex = Regex::new(r"^[!-)+-<>-~][!-~]*$").unwrap();
    let mut end_kmers: Vec<(Kmer, KmerData)> = Vec::new();
    let mut start_kmers: HashMap<Kmer, SmallVec<[KmerData; 4]>> = HashMap::new();
    let reader = Reader::new(input);
    let k = k as usize;
    let mut buffer: Vec<u8> = vec![0; k];
    // None if the record is skipped due to being too short
    let mut identifiers: Vec<Option<String>> = Vec::new();
    for (record_index, record) in reader.records().enumerate() {
        let record = record.context("Failed to parse record from FASTA file")?;
        let seq = record.seq();
        if seq.len() < k.max(min_contig_length) {
            identifiers.push(None);
            continue;
        } else {
            let id = record.id().to_owned();
            if !id_regex.is_match(&id) {
                bail!("Invalid record identifier: {}.\nIdentifier names are restricted by the GFA format.", id);
            }
            identifiers.push(Some(id));
        }
        let (start_fw, end_rv) = if let Some((a, b)) = Kmer::make_pair(&seq[..k], &mut buffer) {
            (a, b)
        } else {
            continue;
        };
        let (end_fw, start_rv) =
            if let Some((a, b)) = Kmer::make_pair(&seq[seq.len() - k..], &mut buffer) {
                (a, b)
            } else {
                continue;
            };
        let fw_data = KmerData::try_new(record_index, false)?;
        let rv_data = KmerData::try_new(record_index, true)?;
        end_kmers.push((end_fw, fw_data));
        end_kmers.push((end_rv, rv_data));
        start_kmers.entry(start_fw).or_default().push(fw_data);
        start_kmers.entry(start_rv).or_default().push(rv_data);
    }
    let mut edges: Vec<Edge> = Vec::new();
    for (end_kmer, end_data) in end_kmers.iter() {
        if let Some(v) = start_kmers.get(end_kmer) {
            for start_data in v.iter() {
                edges.push(Edge {
                    from: *end_data,
                    to: *start_data,
                })
            }
        }
    }
    // Remove unused names to free up memory
    drop(start_kmers);
    drop(end_kmers);
    let mut used = vec![false; identifiers.len()];
    for edge in edges.iter() {
        used[edge.from.index()] = true;
        used[edge.to.index()] = true;
    }
    for (is_used, s) in used.iter().zip(identifiers.iter_mut()) {
        if !is_used {
            *s = None;
        }
    }
    let n_used_elements = identifiers
        .iter()
        .enumerate()
        .rev()
        .find(|(_, n)| n.is_some())
        .map(|(i, _)| i + 1)
        .unwrap_or(0);

    identifiers.truncate(n_used_elements);
    identifiers.shrink_to_fit();
    Ok((identifiers, edges))
}

fn write_gfa(identifiers: &[Option<String>], edges: &[Edge]) -> Result<()> {
    let mut out = BufWriter::new(stdout().lock());
    out.write_all(b"H\tVN:Z:1.2\n")?;
    for edge in edges.iter() {
        out.write_all(b"L\t")?;
        out.write_all(identifiers[edge.from.index()].as_ref().unwrap().as_bytes())?;
        out.write_all(b"\t")?;
        out.write_all(std::slice::from_mut(
            &mut b"+-"[edge.from.is_rc() as usize].clone(),
        ))?;
        out.write_all(b"\t")?;
        out.write_all(identifiers[edge.to.index()].as_ref().unwrap().as_bytes())?;
        out.write_all(b"\t")?;
        out.write_all(std::slice::from_mut(
            &mut b"+-"[edge.to.is_rc() as usize].clone(),
        ))?;
        out.write_all(b"\t*\n")?;
    }
    Ok(())
}

const LONG_ABOUT: &str = "Print a minimal GFA v1.2 file to stdout from a MEGAHIT contig file.
Output file only contains the H (header) line and minimal L (link) lines.
See more information in the README.md.
Usage: megagfa -i final.contigs.fa -k 141 > links.gfa";

#[derive(Parser)]
#[command(author, version, about, long_about = LONG_ABOUT)]
struct Cli {
    /// Value of --k-max used in assembly
    #[arg(short)]
    k: u8,

    /// Input file [stdin if not passed]
    #[arg(short)]
    i: Option<PathBuf>,

    /// Minimum contig length
    #[arg(short, default_value_t = 200)]
    min_contig_length: u32,
}
