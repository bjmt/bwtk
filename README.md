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
GBs of memory for larger files). bwtk is otherwise a fairly lightweight program.

## Installation

```sh
make libBigWig
make release
```

## Quick start

```
bwtk v1.1.0  Copyright (C) 2025  Benjamin Jean-Marie Tremblay
Usage:  bwtk <subcommand> [options]
Available subcommands:
    bw2bg      Convert a bigWig file to bedGraph
    bg2bw      Convert a bedGraph file to bigWig
    merge      Average multiple bigWig files together
    values     Return bigWig values from overlapping BED ranges
    subset     Subset a bigWig using a BED file
    chroms     Print a chrom.sizes file from a bigWig header
    adjust     Perform an operation on bigWig values
    help       Print this message and exit
    version    Print the version number and exit
For subcommand usage, try: bwtk <subcommand> -h
```


