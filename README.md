# nsearch 

NGS data processing library for the rest of us. Written in C++11.

[![CircleCI](https://circleci.com/gh/stevschmid/nsearch.svg?style=shield)](https://circleci.com/gh/stevschmid/nsearch)

## Features

### Database search

- **Query** a list of DNA/RNA/amino acid sequences in a database of your choice.

### Read processing

- **Merge** forward and reverse (Illumina) short-read sequences.
- **Filter** (merged) sequences based on the number of expected errors.

### File formats
Currently supported:

- FASTA (merging output, searching input)
- FASTQ (merging input, merging output, searching input)
- ALNOUT (searching output)

Gzipped input files (e.g. `db.fasta.gz`) are supported.

### Library

Use the `libnsearch` components for your own bioinformatics pipeline.

## Requirements

* cmake: Build system. On macOS, installation via homebrew is recommended: `brew install cmake`.
* (optional) zlib: Reading compressed file support (`.gz`). `brew install zlib`.

## Build instructions

```bash
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make
```

## Usage

Invoke the `nsearch` binary in the `build/nsearch` directory.

## Tests

Run the test in the `build` directory:

```bash
make check
```

## Code Style

A `.clang-format` for [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html) file is provided in the repository. 
