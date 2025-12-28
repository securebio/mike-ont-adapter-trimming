# ont-trim

High-performance ONT adapter trimming tool with streaming and parallelization.

## Features

- **Streaming architecture**: Processes reads in chunks with bounded memory usage
- **Parallel processing**: Near-linear scaling with CPU cores using Rayon
- **Porechop-compatible**: Implements Porechop ABI's alignment and trimming logic
- **Automatic reverse complement**: Searches both orientations of each adapter
- **Read splitting**: Detects and splits chimeric reads with internal adapters
- **Flexible I/O**: Supports stdin/stdout, auto-detects gzip compression

## Installation

```bash
cargo build --release
./target/release/ont-trim --help
```

## Usage

### Basic adapter trimming

```bash
ont-trim -i reads.fastq.gz -o trimmed.fastq.gz \
    --adapter AATGTACTTCGTTCAGTTACGTATTGCT
```

### With multiple adapters

```bash
ont-trim -i reads.fastq.gz -o trimmed.fastq.gz \
    --adapter AATGTACTTCGTTCAGTTACGTATTGCT \
    --adapter GCAATACGTAACTGAACGAAGT
```

### Custom thresholds

```bash
ont-trim -i reads.fastq.gz -o trimmed.fastq.gz \
    --adapter AATGTACTTCGTTCAGTTACGTATTGCT \
    --end-threshold 0.85 \
    --middle-threshold 0.95 \
    -j 8
```

### Streaming from stdin

```bash
cat reads.fastq | ont-trim --adapter ACGT... | gzip > trimmed.fastq.gz
```

## Parameters

### Input/Output
- `-i, --input <FILE>`: Input FASTQ/FASTA file (plain or gzipped)
- `-o, --output <FILE>`: Output file
- `-j, --threads <N>`: Number of threads (default: auto-detect)

### Adapter Detection
- `-a, --adapter <SEQ>`: Adapter sequence (can be specified multiple times)
- `--end-threshold <0.0-1.0>`: Partial identity for end adapters (default: 0.75)
- `--middle-threshold <0.0-1.0>`: Full identity for middle adapters (default: 0.90)

### Trimming Behavior
- `-e, --end-size <BP>`: Region size to search at read ends (default: 150)
- `--min-trim-size <BP>`: Minimum alignment length to trigger trim (default: 4)
- `-m, --min-length <BP>`: Minimum output fragment size (default: 1000)
- `--extra-end-trim <BP>`: Extra bases to trim at ends (default: 2)
- `--extra-middle-good <BP>`: Extra trim on expected side of middle adapter (default: 10)
- `--extra-middle-bad <BP>`: Extra trim on unexpected side (default: 100)
- `--split <true|false>`: Split reads at internal adapters (default: true)

### Performance
- `--chunk-size <N>`: Reads per processing chunk (default: 1000)
- `-c, --compression-level <0-9>`: Gzip compression level (auto-detected from extension)

## Algorithm

### Alignment
Uses **local alignment with affine gap scoring** from the `bio` crate:
- Match: +3
- Mismatch: -6
- Gap open: -5
- Gap extend: -2

This matches Porechop ABI's scoring scheme and handles ONT's "bursty" indel errors effectively.

### Two-Tier Identity Scoring

| Location | Score Type | Threshold | Purpose |
|----------|------------|-----------|---------|
| End adapters | Partial identity | 75% | More lenient for partial adapters at read boundaries |
| Middle adapters | Full identity | 90% | Stricter requirement for chimera detection |

### Extra Trimming

Beyond the detected adapter boundaries, additional bases are removed:
- **End adapters**: +2 bp
- **Middle adapters (good side)**: +10 bp
- **Middle adapters (bad side)**: +100 bp

### Multiple Internal Adapters

Uses iterative masking (like Porechop ABI):
1. Find best adapter match in middle region
2. Mask found adapter
3. Search again until no more matches
4. Split read at all found positions

## Performance

**Note**: These are expected characteristics based on the reference architecture (nao-rustmasker). Actual performance has not yet been benchmarked with real ONT data.

- **Memory**: Expected ~20MB per chunk (1000 reads × ~10kb each) - bounded regardless of file size
- **Parallelization**: Expected near-linear scaling up to 8-16 cores (based on nao-rustmasker's measured performance)
- **Throughput**: I/O (FASTQ parsing, gzip) or alignment could be the bottleneck depending on data characteristics
  - The `bio` crate uses scalar (non-SIMD) alignment, which may be slower than expected
  - If profiling shows alignment is the bottleneck, can swap to `parasail-rs` for SIMD acceleration

## Implementation Notes

### Architecture Pattern (from nao-rustmasker)
```
Read → Chunk (1000 reads) → Parallel Process → Sequential Write → Repeat
```

### Key Dependencies
- `needletail`: FASTQ/FASTA parsing with auto-gzip detection
- `bio`: Sequence alignment with affine gap scoring
- `rayon`: Data parallelism
- `gzp`: Parallel gzip compression
- `clap`: CLI argument parsing

## Future Enhancements

- Adapter inference (learn from data like Porechop_ABI)
- Quality-aware trimming boundaries
- Statistics output (adapter frequency, trim distribution)
- JSON/TSV report mode
- Barcode handling

## References

- [Porechop ABI](https://github.com/bonsai-team/Porechop_ABI) - Original Python implementation
- [nao-rustmasker](../nao-rustmasker) - Reference architecture for streaming/parallelization
