# megagfa
Create a minimal [GFA file](https://github.com/GFA-spec/GFA-spec) from an intermediate MEGAHIT contig file, with links between contigs that are linked in the assembly graph.

## Background
The metagenomic assembler [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) is based on relative simple techniques.
In particular, it's based on the original iterative de Bruijn graph approach with few modifications.
In this paradigm, contigs are equal every maximally-long path in the de Bruijn graph where every vertex is visited at most once, and where every vertex except possibly the last has outdegree one.
Since every vertex is a kmer, this implies that if contig A is a continuation of contig B, then the first kmer is contig A is the last kmer of contig B.

In other words, if contig A ends with the same kmer contig B begins with, there is a high probability that contig B is the true continuation of contig A. Though of course not high enough probability the assembler is willing to fuse the two contigs.
This information is useful for binning.

## Installation
* [Install Rust](https://www.rust-lang.org/tools/install)
* Navigate to this directory, then compile with `cargo build --release`.
* The binary can be found in `target/release/megagfa`

To compile for Computerome, use `cargo build --release --target=x86_64-unknown-linux-musl`

## How to use
```
$ megagfa -i final.contigs.fa -k 141 > links.gfa
```
* The value passed with `-k` must be equal to the value of `--k-max` that MEGAHIT was run with, else the results will be wrong.
  The default value for MEGAHIT is 141.
  You can find the value for any given run as the largest number in the "k list" printed in the log file.
* If `-i` is not passed, the program will read from stdin. Hence, to read compressed contigs, use `gunzip -dc final.contigs.fa.gz | megagfa -k 141 > links.gfa`.

## Output
The output printed to stdout is a GFA 1.2 file. It looks like this:
```
H       VN:Z:1.2
L       k141_100502     +       k141_11333      -       *
L       k141_100502     +       k141_31468      -       *
L       k141_0  +       k141_110078     -       *
L       k141_0  +       k141_48545      +       *
L       k141_33502      +       k141_1046       -       *
```
First it outputs the header saying it's a GFA version so-and-so.
Then, each line beginning with `L` specifies a shared k-mer that constitute an edge in the assembly graph between two contigs. The plus and minus means forward / reverse strand, respectively.
Hence, the first `L` line says that contig `k141_100502` ends with the same k-kmer that the reverse-complement of `k141_11333` starts with.
