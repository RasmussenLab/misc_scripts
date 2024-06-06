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
    // We can read from stdin, from a file, or from a gzipped file. In any case, we wrap the result
    // in a BufRead so we can guarantee the input type implements BufRead.
    let input: Box<dyn BufRead> = if let Some(p) = args.i {
        if !p.is_file() {
            exitwith(&format!(
                "Not an existing file: \"{}\"",
                p.to_string_lossy()
            ));
        }
        // Check if the user passes e.g. a file k79.contigs.fa, but passes -k 75, where the values
        // of k differ. This will raise an error.
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
        // Return a BufReader wrapping either the opened file directly, or a gzip reader if the file name
        // ends with .gz.
        let file = std::fs::File::open(p.clone())
            .with_context(|| format!("Could not open input file \"{}\"", p.to_string_lossy()))?;
        if p.extension().is_some_and(|e| e == "gz") {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    } else {
        // A locked stdin implements BufRead, so we can use that directly.
        // Also, no other part of this program reads stdin, so there is no downside.
        Box::new(stdin().lock())
    };
    let (identifiers, edges) = find_edges(input, args.k, args.min_contig_length as usize)?;
    print_gfa(&identifiers, &edges)?;
    Ok(())
}

// Contains the index of the contig the kmer came from, and whether it's reverse-complement
// or not.
// The last bit of information needed to identify a kmer is whether it's the ending or the starting kmer,
// but this is stored implicitly as the end/start kmers are stored in two different
// data structures in this program.
// The information is packed into 32 bits in order to save memory, and to make the data structures
// that store KmerData smaller and therefore faster.
#[derive(Clone, Copy)]
struct KmerData(u32);

impl KmerData {
    fn try_new(index: usize, rc: bool) -> Result<Self> {
        // It's unlikely we get more than 2 billion records in a file, but let's check it anyway
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

// From: The ending kmer. To: The starting kmer of the next contig.
struct Edge {
    from: KmerData,
    to: KmerData,
}

// We use this LUT (lookup table) to encode arbitrary DNA/RNA nucleotides into two bits.
// This is to make the Kmer struct smaller - both for memory reasons, but also to
// make hashing it faster.
// This will mean contigs with ambiguous nucleotides in the start/ending kmers will be skipped,
// but I'm not sure MEGAHIT can even process ambiguous kmers in its graph anyway, so no loss.
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

// See comments above make_lut
#[derive(PartialEq, Eq, Hash)]
struct Kmer(Vec<u8>);

impl Kmer {
    // Return the (forward, reverse_complement) pair, or None if the kmer
    // contains any bytes not in "AaCcGgTtUu".
    // The into is a temporary buffer that can be reused between calls.
    fn make_pair(bytes: &[u8], into: &mut [u8]) -> Option<(Self, Self)> {
        translate(bytes, into)?;
        let fw = compact(into);
        reverse_complement(into);
        let rv = compact(into);
        Some((Kmer(fw), Kmer(rv)))
    }
}

// Encode a slice of input DNA/RNA bytes to a slice of encodings, where the encodings
// are the two lower bits of the byte.
fn translate(s: &[u8], into: &mut [u8]) -> Option<()> {
    let _: [u8; 256] = LUT; // Statically check LUT has 256 elements
    for (byte, to) in s.iter().zip(into.iter_mut()) {
        // Safety: a 256-length array cannot be indexed out of bounds with a u8.
        unsafe {
            let b = *LUT.get_unchecked(*byte as usize);
            if b == 0xff {
                return None;
            } else {
                *to = b
            }
        }
    }
    Some(())
}

// This works if x is encoded using `translate`.
// XORing with 3 flips the last 2 bits, switching the encoding of A<->T and C<->G.
fn reverse_complement(x: &mut [u8]) {
    x.reverse();
    for i in x.iter_mut() {
        *i ^= 3
    }
}

// From a slice of encodings where only the two lowest bits are used,
// we can compact it to a 4x smaller vector.
// TODO: This ought to be better using as_chunks, but this is not stabilized
// as of writing.
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

// Compact a slice of 0-4 encodings to a single u8.
fn compact_element(s: &[u8]) -> u8 {
    s.iter().fold(0, |old, new| (old << 2) | new)
}

fn find_edges(
    input: impl BufRead,
    k: u8,
    min_contig_length: usize,
) -> Result<(Vec<Option<String>>, Vec<Edge>)> {
    // According to the GFA specs, FASTA identifiers must conform to this regex.
    // Too bad if we have identifiers which don't - we must end the program.
    let id_regex = Regex::new(r"^[!-)+-<>-~][!-~]*$").unwrap();

    // Approach: We store the end kmers in a simple vector, then all the starting kmers in a
    // hash map hashed by the kmer itself.
    // After we've collected start/end kmers, we loop over every end kmer, find all the start kmers
    // that have the same Kmer content, then create an edge from end -> start.
    let mut end_kmers: Vec<(Kmer, KmerData)> = Vec::new();
    let mut start_kmers: HashMap<Kmer, SmallVec<[KmerData; 2]>> = HashMap::new();
    let reader = Reader::new(input);
    let k = k as usize;
    // Temporary buffer used for calls to make_pair.
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
        // The start of the sequence contains the starting forward kmer and the ending reverse kmer.
        // The continue here happens if the seq contains symbols that cannot be enocded to a kmer, which
        // I believe is rare an relatively unimportant to handle.
        let Some((start_fw, end_rv)) = Kmer::make_pair(&seq[..k], &mut buffer) else {
            continue;
        };

        // Vice versa - the end of the sequence contains the ending forward kmer and the starting reverse kmer
        let Some((end_fw, start_rv)) = Kmer::make_pair(&seq[seq.len() - k..], &mut buffer) else {
            continue;
        };
        let fw_data = KmerData::try_new(record_index, false)?;
        let rv_data = KmerData::try_new(record_index, true)?;
        end_kmers.push((end_fw, fw_data));
        end_kmers.push((end_rv, rv_data));
        start_kmers.entry(start_fw).or_default().push(fw_data);
        start_kmers.entry(start_rv).or_default().push(rv_data);
    }
    // Now, for every end kmer, we see if there are any matching starting kmers, then
    // we create an edge from end kmer to start kmer.
    // Why not from start to end? Remember, if contig B follows contig A, then we
    // go from the last contig of A to the first contig of B.
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
    // Free up unneeded memory - we can delete all identifiers that have no edge,
    // then shrink the identifiers vector to not contain trailing Nones.
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

fn rc_byte(rc: bool) -> &'static [u8] {
    if rc {
        b"-"
    } else {
        b"+"
    }
}

// Write a minimal GFA
fn print_gfa(identifiers: &[Option<String>], edges: &[Edge]) -> Result<()> {
    let mut out = BufWriter::new(stdout().lock());
    // Write header - this is GFA version 1.2
    out.write_all(b"H\tVN:Z:1.2\n")?;
    for edge in edges.iter() {
        // Write L lines: L
        out.write_all(b"L\t")?;
        // Name of sequende with end kmer (from)
        out.write_all(identifiers[edge.from.index()].as_ref().unwrap().as_bytes())?;
        out.write_all(b"\t")?;
        // Whether the from sequence is forward or reverse
        out.write_all(rc_byte(edge.from.is_rc()))?;
        out.write_all(b"\t")?;
        // Same for the to edge
        out.write_all(identifiers[edge.to.index()].as_ref().unwrap().as_bytes())?;
        out.write_all(b"\t")?;
        out.write_all(rc_byte(edge.to.is_rc()))?;
        // A star for the missing overlap (which carries no information, the user should know
        // it's always just one kmer's overlap)
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
