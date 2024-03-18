# dnazip
This program walks through a given directory, gzip-compressing any uncompressed FASTA or FASTQ files.

Examples:
* `dnazip --threads 15 --verbose my_dir`
* `dnazip --dry-run .`

Does not follow symbolic links (so no infinite loops)
