# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TagDust2 is a bioinformatics tool for processing next-generation sequencing (NGS) data. It extracts and labels sequences containing adapter, linker, barcode, and fingerprint sequences using Hidden Markov Models (HMMs). The tool can demultiplex reads, automatically detect library preparation methods, and filter out contaminants.

## Build System

This project uses GNU Autotools for configuration and building:

### Essential Build Commands
```bash
# Initial setup (run once after cloning)
./autogen.sh
./configure

# Build the project
make

# Run tests
make check

# Install system-wide
make install
```

### Development Build Options
```bash
# Enable debugging
./configure --enable-debugging

# Enable valgrind tests
./configure --enable-valgrind-tests

# Build with specific compiler flags
./configure CFLAGS="-O2 -Wall -std=gnu99"
```

### Testing
- Main test suite: `make check`
- Individual test programs are built in `src/`: `tagdustiotest`, `tagdust_rtest`, `simreads_rtest`, `evalres_rtest`
- Development tests in `dev/`: `sanity_test.sh`, `casava_test.sh`, `bar_read_test.sh`

## Architecture

### Core Components
- **main.c**: Main tagdust program entry point
- **barcode_hmm.c/h**: HMM construction, training, and sequence matching
- **interface.c/h**: Command-line argument parsing and program initialization
- **io.c/h**: File I/O operations for FASTQ/FASTA formats
- **nuc_code.c/h**: Nucleotide encoding and sequence operations
- **misc.c/h**: Utility functions and data structures
- **kslib.c/h**: Sequence parsing and utility library

### Executables Built
1. **tagdust**: Main demultiplexing tool
2. **simreads**: Read simulation utility
3. **evalres**: Results evaluation tool
4. **merge**: Sequence merging utility
5. **rename_qiime**: QIIME format renaming tool

### Key Data Structures
- HMM structures for sequence architecture modeling
- Sequence data structures with quality scores
- Barcode matching and scoring systems

## Configuration Files

### Architecture Files
Architecture files define the expected structure of reads:
- Located in `dev/` and `casava_demo/` directories
- Examples: `casava_arch.txt`, various `EDITTAG_*` files
- Format: Specify barcodes, read segments, and expected patterns

### Test Data
- `dev/`: Contains test FASTQ files and expected output files
- `benchmark/`: Performance testing scripts and data
- Gold standard files ending in `_gold.txt` for validation

## Development Workflow

### Adding New Features
1. Modify relevant source files in `src/`
2. Update `src/Makefile.am` if adding new source files
3. Run `make` to build
4. Run `make check` to ensure tests pass
5. Add new tests in `dev/` if needed

### Debugging
- Configure with `--enable-debugging` for debug builds
- Use `--enable-debugging=1|2|3` for different verbosity levels
- Debug macros available via `debug_print()` in `tagdust2.h`

### Memory Management
- Custom memory allocation macros in `malloc_macro.h`
- SIMD-aligned memory allocation support
- Memory debugging available with debug builds

## Paper Reproducibility

To reproduce the benchmarks from the TagDust2 paper:
```bash
cd reproducibility/scripts
make -f run.mk benchmark
```

This requires 64-bit system and may take considerable time.