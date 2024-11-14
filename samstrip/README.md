# samstrip
This program removes parts of SAM files that do not give alignment information, saving disk space.
The file is read from stdin and printed to stdout.
The SEQ and QUAL fields are zeroed, and all optional fields (except the NM field, which is not really opional) is removed.

Usage:
The following code converts a BAM file to a SAM file, strips it, then converts it back.
`samtools view -h file.bam | samstrip | samtools view -b > stripped.bam`
