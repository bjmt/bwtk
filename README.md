# bwtk: A bigWig Toolkit

bwtk is a program containing a set of utilities for handling [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html)
files, as well as converting to and from the bedGraph format. bigWig files
are indexed binary files which allow for fast read access of different
parts of the file, ideal for genome browsers and extracting data from specific
chromosomes or subsequences quickly. It has some limitations however, such as
only being able to work with genomes that don't contain any individual
chromosome larger than the size limit of `uint32_t` (around 4.2 billion).
This program is built upon the
[libBigWig](https://github.com/dpryan79/libBigWig) library to read and write
bigWigs. One note of caution regarding bwtk's usage of this library: during the
final stages of creating new bigWigs, an indexing step occurs which can require
a very large amount of memory (on the order of hundreds of MBs, perhaps even
GBs of memory for larger files). Therefore, expect any subcommand which involves creating new bigWig files to take up a decent chunk of memory. bwtk is otherwise a fairly lightweight and fast program, needing only a few MBs of memory for the remaining subcommands.

bwtk mostly reimplements some of the existing functionality from the original UCSC tools, with a few additions and improvements:

- bedGraph files can be gzipped when converting to bigWig
- Multiple bigWigs can be averaged/summed/min'd/max'd together
- Retrieval of single-base resolution data from BED ranges
- Subsetting of bigWigs
- Value operations such as addition, multiplication, log10 transformation
- Binning of bigWig values and collapsing of any resulting sequential ranges with identical values, leading to substantially smaller file sizes depending on the bin step value

## Installation

```sh
make libz libBigWig
make release
```

## Quick start

```
bwtk v1.5.0  Copyright (C) 2025  Benjamin Jean-Marie Tremblay
Usage:  bwtk <subcommand> [options]
Available subcommands:
    bg2bw      Convert a bedGraph file to bigWig
    adjust     Perform an operation on a bigWig
    merge      Average multiple bigWig files together
    values     Return bigWig values from overlapping BED ranges
    score      Get summary scores of bigWig values from BED ranges
    chroms     Print a chrom.sizes file from a bigWig header
    help       Print this message and exit
    version    Print the version number and exit
For subcommand usage, try: bwtk <subcommand> -h
```

## File transformations

- `bg2bw`: one bedGraph -> [optional operation] -> one bigWig
- `adjust`: one bigWig -> [optional operation] -> one bigWig (or bedGraph)
- `merge`: Multiple bigWigs -> [merge values] -> [optional operation] -> one bigWig

## Data extraction

- `values`: Extract single base-resolution range scores
- `score`: Calculate summary information of range scores
- `chroms`: Retrieve chromosome names and sizes


