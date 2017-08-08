# nsearch (WIP title)

Metagenomics tool/library for the rest of us. Written in C++11.

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

### Library

Use the `libnsearch` components for your own bioinformatics pipeline.

## Requirements

You need `cmake` (build system). On macOS, installation via homebrew is recommended:

```bash
brew install cmake
```

## Build instructions

```bash
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make
```

## Usage

Invoke the `app/nsearch` binary in the `build` directory.

## Tests

Run the test in the `build` directory:

```bash
make check
```

## Code Style

A `.clang-format` for [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html) file is provided in the repository. 
