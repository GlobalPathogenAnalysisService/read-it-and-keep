# read-it-and-keep
Read contamination removal.


## Installation
### From source
Make the executable `src/readItAndKeep` by running:
```
cd src && make
```

### Bioconda ![Platforms](https://anaconda.org/bioconda/read-it-and-keep/badges/platforms.svg)

From an existing environment:
```
conda install -c bioconda read-it-and-keep
```
Using a new environment (recommended):
```
conda create -n read-it-and-keep -c bioconda python=3 read-it-and-keep
conda activate read-it-and-keep
```

### Docker
Build a docker container by cloning this repository and running:
```
docker build -f Dockerfile -t <TAG> .
```

### Singularity
Build a singularity container by cloning this repository and running:
```
sudo singularity build readItAndKeep.sif Singularity.def
```


## Usage

To run on paired Illumina reads, in two files `reads1.fq.gz` and `reads2.fq.gz`:

```
readItAndKeep --ref_fasta ref_genome.fasta --reads1 reads1.fq.gz --reads2 reads2.fq.gz -o out
```

This will output `out.reads_1.fastq.gz` and `out.reads_2.fastq.gz`.

To run on one file of nanopore reads `reads.fq.gz`:

```
readItAndKeep --tech ont --ref_fasta ref_genome.fasta --reads1 reads.fq.gz -o out
```

This will output `out.reads.fastq.gz`.

If the input reads files are in FASTA format, then it will output reads in FASTA format, calling the files `*.fasta.*` instead of `*.fastq.*`.

It always writes the counts of input and output reads to `STDOUT` in tab-delimited format, for example:

```
Input reads file 1	1000
Input reads file 2	1000
Kept reads 1	950
Kept reads 2	950
```

All logging messages are sent to `STDERR`.

**Required arguments:**

- `--ref_fasta`: reference genome in FASTA format.
- `--reads1`: at least one reads file in FASTA[.GZ] or FASTQ[.GZ] format.
- `-o,--outprefix`: prefix of output files.

Please note there is an option `--tech`, which defaults to `illumina`. Use `--tech ont` for nanopore reads.

**Optional arguments:**

- `--reads2`: name of second reads file, i.e. mates file for paired reads
- `--enumerate_names`: rename the reads `1`,`2`,`3`,... (for paired reads, will also add `/1` or `/2` on the end of names)
- `--debug`: debug mode. More verbose and writes debugging files
- `--min_map_length`: minimum length of match required to keep a read in bp (default `50`)
- `--min_map_length_pc`: minimum length of match required to keep a read, as a percent of the read length (default `50.0`)
- `-V,--version`: show version and exit

### Docker
Additional arguments need to be supplied to allow Docker to access input and output files. Below is a functional example:

```
docker run /path/to/read-it-and-keep/tests:/tests  [-v /path/to/input:/input -v /path/to/output:/output] <TAG> --ref_fasta /tests/MN908947.3.fa --reads1 /input/<SAMPLE>_1.fastq.gz --reads2 /input/<SAMPLE>_2.fastq.gz --outprefix /output/
```
## Tests

These are under development. To run them you will need:
1. Python 3
2. Python package [pytest](https://docs.pytest.org/en/stable/) (`pip install pytest`)
3. Python package [pyfastaq](https://github.com/sanger-pathogens/Fastaq)  (`pip install pyfastaq`)
4. [ART read simulator](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
   installed, so that `art_illumina` is in your `$PATH`
5. [badread](https://github.com/rrwick/Badread) for nanopore read simulation.

Run the tests after compiling the source code, ie:
```
cd src
make
make test
```

## Acknowledgements

This repository includes unedited copies of the code from:
* [gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/), [LGPL 2.1 licence](https://github.com/GenomePathogenAnalysisService/read-it-and-keep/blob/main/src/ext/gzstream/COPYING.LIB)
* [minimap2](https://github.com/lh3/minimap2), [MIT licence](https://github.com/GenomePathogenAnalysisService/read-it-and-keep/blob/main/src/ext/minimap2-2.22/LICENSE.txt)
* [CLI11](https://github.com/CLIUtils/CLI11) header file, licence is at start of [cli11.hpp](https://github.com/GenomePathogenAnalysisService/read-it-and-keep/blob/main/src/CLI11.hpp)
