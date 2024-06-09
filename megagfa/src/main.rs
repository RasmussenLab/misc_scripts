use anyhow::{self, bail, Context, Result};
use bio::io::fasta::Reader;
use clap::Parser;
use smallvec::SmallVec;
use std::{
    collections::HashMap,
    io::{stdin, stdout, BufRead, BufReader, BufWriter, Write},
    num::NonZeroU8,
    path::PathBuf,
};

fn exitwith(s: &str) -> ! {
    eprintln!("{}", s);
    std::process::exit(1)
}

fn main() -> Result<()> {
    let args = Cli::parse();

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
            s.strip_prefix('k').and_then(|s| s.strip_suffix(".contigs.fa")).and_then(|s| {
                s.parse::<u8>().ok().map(|k| {
                    if k != args.k.get() {
                        exitwith(&format!("ERROR: K value passed with -k is {}, but given file is {} with different K value.", args.k, k))
                    };
                })
            })
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
// that store KmerOrigin smaller and therefore faster.
#[derive(Clone, Copy, Debug)]
struct KmerOrigin(u32);

impl KmerOrigin {
    fn is_rc(&self) -> bool {
        self.0 > 0x7fffffff
    }

    fn index(&self) -> usize {
        (self.0 & 0x7fffffff) as usize
    }

    fn reverse_complement(&self) -> Self {
        Self(self.0 ^ 0x80000000)
    }
}

// Just a convenience struct so we can construct a fw and an rc KmerOrigin in one go
struct KmerOriginPair {
    fw: KmerOrigin,
    rc: KmerOrigin,
}

impl KmerOriginPair {
    fn try_new(index: usize) -> Result<Self> {
        // It's unlikely we get more than 2 billion records in a file, but let's check it anyway
        let x: u32 = index
            .try_into()
            .ok()
            .and_then(|u| if u > 0x7fffffff { None } else { Some(u) })
            .context("Can only hande 2147483647 FASTA records in one file")?;
        Ok(Self {
            fw: KmerOrigin(x),
            rc: KmerOrigin(x | 0x80000000),
        })
    }
}

// From: The ending kmer. To: The starting kmer of the next contig.
struct Edge {
    from_end: KmerOrigin,
    to_start: KmerOrigin,
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

// We encode the kmers in two bits, so this is identical to ceiling dividing by 4.
fn encoding_size(k: NonZeroU8) -> NonZeroU8 {
    // Safety: If k > 3, the first term is > 0 and < 254.
    // If k is 1, 2, or 3, the second term is > 0.
    // Hence the sum will always be nonzero and cannot overflow.
    unsafe {
        (k.get() / 4 + ((k.get() % 4) > 0) as u8)
            .try_into()
            .unwrap_unchecked()
    }
}

// Dense representation of all observed kmers, packed into a single vector.
struct Kmers {
    mers: Vec<u8>, // The DNA kmers themselves, packed together. Length: k * data.len()
    data: Vec<KmerOrigin>,
    encoding_buffer: Vec<u8>, // length encoding_size(k). Ephemeral.
    k: NonZeroU8,
}

impl Kmers {
    fn new(k: NonZeroU8) -> Self {
        // Preallocate to avoid unnecessary reallocations. This costs about 3 MB up front.
        let assumed_kmers = 200_000;
        Self {
            mers: Vec::with_capacity(encoding_size(k).get() as usize * assumed_kmers),
            data: Vec::with_capacity(assumed_kmers),
            encoding_buffer: vec![0; encoding_size(k).get() as usize],
            k,
        }
    }

    // How to get the kmers and KmerOrigin out of this struct.
    fn iter_kmers(&self) -> impl Iterator<Item = (&KmerOrigin, &[u8])> {
        let chunk_size = self.mers.len() / self.data.len();
        self.data.iter().zip(self.mers.chunks_exact(chunk_size))
    }

    // Add the kmers and kmer data from a sequence to this struct.
    // None if seq too small, or contains non-DNA
    fn add(&mut self, seq: &[u8], index: usize) -> Option<()> {
        let KmerOriginPair {
            fw: fwdata,
            rc: rvdata,
        } = KmerOriginPair::try_new(index).ok()?;

        // Add forward starting kmer
        if translate(
            seq.get(0..self.k.get() as usize)?,
            &mut self.encoding_buffer,
        )
        .is_some()
        {
            self.mers.extend_from_slice(&self.encoding_buffer);
            self.data.push(fwdata);
        }

        // Add reverse starting kmer
        translate(
            seq.get(seq.len() - self.k.get() as usize..)?,
            &mut self.encoding_buffer,
        )?;
        reverse_complement(self.k, &mut self.encoding_buffer);
        self.mers.extend_from_slice(&self.encoding_buffer);
        self.data.push(rvdata);
        Some(())
    }
}

// None if the sequence contains a byte which are not ACGTUacgtu.
fn translate(seq: &[u8], into: &mut [u8]) -> Option<()> {
    // Handle first chunks of 4, which each are translated to a single byte.
    let chunks = seq.chunks_exact(4);
    let (last_encoding, mut is_error) = translate_chunk(chunks.remainder().iter());
    for (e, (b, err)) in into
        .iter_mut()
        .zip(chunks.map(|chunk| translate_chunk(chunk.iter())))
    {
        *e = b;
        is_error |= err;
    }
    // Handle last element. We could handle all elements in a single loop using
    // chunk instead of chunks_exact, but that would cause worse code to be emitted.
    if let Some(e) = into.last_mut() {
        *e = last_encoding
    };
    if is_error {
        None
    } else {
        Some(())
    }
}

fn translate_chunk<'a, T: Iterator<Item = &'a u8>>(x: T) -> (u8, bool) {
    // Statically and locally verify it has 256 elements for safety
    let lut: [u8; 256] = LUT;
    x.fold((0, false), |(kmer, is_error), byte| unsafe {
        let b = lut.get_unchecked(*byte as usize);
        ((kmer << 2) | b, is_error | (*b == 0xff))
    })
}

fn reverse_complement(k: NonZeroU8, kmer: &mut [u8]) -> &[u8] {
    // First we reverse. We need to reverse each byte (chunk of 4 2-bit symbols)
    // then we bitreverse each byte.
    // So e.g. a byte like ABCDEFGH becomes HGFEDCBA, when it should be
    // GHEFCDAB. So, we also do a bit of bitshuffling to get correctly reversed.
    // It should optimise well.
    kmer.reverse();
    for i in kmer.iter_mut() {
        *i = {
            let r = i.reverse_bits();
            let bitreversed = ((r & 0b10101010) >> 1) | ((r & 0b01010101) << 1);
            // Also complement by bitwise not on the bits. This works due to how
            // the nucleotides are stored
            !bitreversed
        }
    }

    // The reversing operation have also reversed where the unused padding bits
    // are. E.g. for a 5-mer it's encoded as AABBCCDD xxxxxxEE, then when reversed
    // its EExxxxxx DDCCBBAA, when the correct result is EEDDCCBB xxxxxxAA.
    // We solve this by shifting the bits.
    let used_bits = 2 * (k.get() % 4);
    let unused_bits = 8 - used_bits;
    let fst = kmer.first_mut().unwrap();
    // First, shift the first byte. In the example above, it's the EExxxxxx shifted by 6.
    *fst >>= unused_bits;
    // It's now xxxxxxEE DDCCBBAA

    // We now need to shift all bytes leftward.
    for i in 0..kmer.len() - 1 {
        // Each byte is its own content shifted leftward, OR'd with the next byte,
        // shifted rightwards.
        unsafe {
            let v =
                (kmer.get_unchecked(i) << unused_bits) | (kmer.get_unchecked(i + 1) >> used_bits);
            *(kmer.get_unchecked_mut(i)) = v;
        }
    }
    // It's now EEDDCCBB DDCCBBAA

    // The last byte just needs its upper bits masked
    unsafe { *kmer.last_mut().unwrap_unchecked() &= (1u8 << used_bits).wrapping_sub(1) };

    // We now have EEDDCCBB xxxxxxAA, the correct answer.
    kmer
}

#[cfg(test)]
mod test_rc {
    use crate::{encoding_size, reverse_complement, translate};
    use std::num::NonZeroU8;

    #[test]
    fn test_rc_fn() {
        for (i, j) in [(b"atcgactacG", b"cGTAGTCGAT")] {
            let n: NonZeroU8 = i
                .len()
                .try_into()
                .ok()
                .and_then(|b| NonZeroU8::new(b))
                .unwrap();
            assert_eq!(n.get() as usize, j.len());
            let mut a = vec![0u8; encoding_size(n).get() as usize];
            let mut b = a.clone();
            translate(i, &mut a).unwrap();
            reverse_complement(NonZeroU8::new(10).unwrap(), &mut a);
            translate(j, &mut b);
            assert_eq!(a, b);
        }
    }
}

// According to the GFA specs, FASTA identifiers must conform to this pattern.
// Too bad if we have identifiers which don't - we must end the program.
fn is_acceptable_identifier(s: &[u8]) -> bool {
    s.split_first().is_some_and(|(first, rest)| {
        ((b'!'..=b')').contains(first)
            | (b'+'..=b'<').contains(first)
            | (b'>'..=b'~').contains(first))
            & rest
                .iter()
                .fold(true, |acc, b| acc & (b'!'..=b'~').contains(b))
    })
}

fn find_edges(
    input: impl BufRead,
    k: NonZeroU8,
    min_contig_length: usize,
) -> Result<(Vec<Option<String>>, Vec<Edge>)> {
    // Approach: We store the starting kmers (forward and reverse-complement)
    // in a HashMap, with keys being kmers and values being KmerOrigin to show
    // where the kmer is from.
    // We can then look up in the hash map to match KmerOrigins with shared kmers
    // and create edges between them
    let mut kmers = Kmers::new(k);
    let reader = Reader::new(input);
    // None if the record is skipped due to being too short
    let mut identifiers: Vec<Option<String>> = Vec::new();
    for (record_index, record) in reader.records().enumerate() {
        let record = record.context("Failed to parse record from FASTA file")?;
        let seq = record.seq();
        if seq.len() >= min_contig_length && kmers.add(record.seq(), record_index).is_some() {
            let id = record.id();
            if !is_acceptable_identifier(id.as_bytes()) {
                bail!("Invalid record identifier: {}.\nIdentifier names are restricted by the GFA format to regex [!-)+-<>-~][!-~]*.", id);
            }
            identifiers.push(Some(id.to_owned()))
        } else {
            identifiers.push(None);
        }
    }
    // Now, for every end kmer, we see if there are any matching starting kmers, then
    // we create an edge from end kmer to start kmer.
    // Why not from start to end? Remember, if contig B follows contig A, then we
    // go from the last contig of A to the first contig of B.
    let mut map: HashMap<&[u8], SmallVec<[KmerOrigin; 2]>> = HashMap::new();
    for (start_data, kmer) in kmers.iter_kmers() {
        map.entry(kmer).or_default().push(*start_data);
    }

    let mut edges: Vec<Edge> = Vec::new();
    // Since `map` borrows from `kmers`, we can't mutate the buffer inside `kmers`
    // and must allocate a new one. No worries.
    let mut rc_buffer: Vec<u8> = vec![0; encoding_size(k).get() as usize];
    for (rc_end_kmer, rc_end_datas) in map.iter() {
        // The map contains starting kmers. By reverse-complementing them, we get
        // ending kmers, which we then use to look up into the map.
        rc_buffer.copy_from_slice(rc_end_kmer);
        let end_kmer = reverse_complement(k, &mut rc_buffer);
        if let Some(start_datas) = map.get(end_kmer) {
            for start_data in start_datas.iter() {
                for rc_end_data in rc_end_datas.iter() {
                    edges.push(Edge {
                        from_end: rc_end_data.reverse_complement(),
                        to_start: *start_data,
                    })
                }
            }
        }
    }
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
        out.write_all(
            identifiers[edge.from_end.index()]
                .as_ref()
                .unwrap()
                .as_bytes(),
        )?;
        out.write_all(b"\t")?;
        // Whether the from sequence is forward or reverse
        out.write_all(rc_byte(edge.from_end.is_rc()))?;
        out.write_all(b"\t")?;
        // Same for the to edge
        out.write_all(
            identifiers[edge.to_start.index()]
                .as_ref()
                .unwrap()
                .as_bytes(),
        )?;
        out.write_all(b"\t")?;
        out.write_all(rc_byte(edge.to_start.is_rc()))?;
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
    k: NonZeroU8,

    /// Input file (may be gzipped) [stdin]
    #[arg(short)]
    i: Option<PathBuf>,

    /// Minimum contig length
    #[arg(short, default_value_t = 200)]
    min_contig_length: u32,
}
