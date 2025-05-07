# SemiBin reclustering
This is an adaptation of the reclustering step in SemiBin. The code is mostly copied from the original SemiBin implementation (which is under an MIT license), but has been adapted to be run from command line, on Vamb data.

## Installation
This script requires a virtual environment, both to make sure certain Python packages are present, and also to ensure the programs `hmmsearch` and `prodigal` can be run from command line.

I recommend installing Pixi (https://pixi.sh) and using that to manage the virtual env.
Use the `pixi.toml` file in this directory to instantiate the Pixi environment.
Even if not using Pixi, you can get a list of dependencies and versions from that toml file.

## Running
You need:

* FASTA file with the contigs
* Vamb / taxvamb etc latent NPZ file
* Vamb / taxvamb clusters TSV file
* For long reads, a TaxVamb taxonomy file

Then:

1. Extract FASTA headers to a file using `grep -E "^>" binbench_cli/data/contigs_2kbp.fna | cut -c 2- > headers.txt`

2. Choose the algorithm: kmeans for short read, dbscan for long reads

3. First time you run it, it will run prodigal and hmmsearch in parallel in the output directory, and produce a `markers.hmmout`.
   This will be super slow. Use lots of threads for this, it scales well. For the following runs, pass this file as `--hmmout_path` to speed up significantly.
   When running hmmsearch, the file `marker.hmm` must be present in the same directory as the `reclustering.py` script.
   This file can be found in the Vamb repo in `vamb/vamb/marker.hmm`.

4. Make the output directory - this script does not make it itself

5. Then run from command line with `pixi run start {ARGUMENTS}`, e.g. `pixi run start vae_clusters_split.tsv latent.npz contigs_2kbp.fna headers.txt outputdir kmeans`.