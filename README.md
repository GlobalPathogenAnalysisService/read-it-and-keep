# read-it-and-keep
Read contamination removal.


## Compile
Make the executable `readItAndKeep`  in the `src` directory by running:
```
cd src && make
```

## Usage
Required options:
1. `--ref_fasta`: reference genome in FASTA format.
2. `--reads1`: at least one reads file in FASTQ[.GZ] format.
3. `-o,--outprefix`: prefix of output files.


Run on one file of reads:
```
readItAndKeep --ref_fasta ref_genome.fasta --reads1 reads1.fq.gz -o out
```
It will output `out.reads_1.fastq.gz`.

Run on two files of reads by adding the option `--reads2 reads2.fq.gz` to the
previous command. It will output `out.reads_1.fastq.gz` and
`out.reads_2.fastq.gz`.

It always writes the counts of input and output reads to `STDOUT` in
tab-delimited format, for example:
```
Input reads file 1	1000
Input reads file 2	1000
Kept reads 1	950
Kept reads 2	950
```
All logging messages sent to `STDERR`.
