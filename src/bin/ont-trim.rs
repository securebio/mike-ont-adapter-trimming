//! ONT adapter trimming CLI tool
//!
//! High-performance streaming adapter trimmer for Oxford Nanopore reads.

use std::fs::File;
use std::io::{BufWriter, Write};

use clap::Parser;
use gzp::{deflate::Gzip, par::compress::ParCompressBuilder, Compression as GzpCompression};
use needletail::{parse_fastx_file, parse_fastx_stdin};
use rayon::prelude::*;

use ont_trim::{
    prepare_adapters, process_read, Adapter, FastqRecord, ProcessedRead, ProcessingConfig,
};

#[derive(Parser, Debug)]
#[command(author, version, about = "High-performance ONT adapter trimming tool")]
struct Args {
    /// Input FASTQ/FASTA file (plain or gzipped). Reads from stdin if not specified.
    #[arg(short = 'i', long)]
    input: Option<String>,

    /// Output file. Writes to stdout if not specified.
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// Adapter sequences to search for (can be specified multiple times).
    /// Both forward and reverse complement are searched automatically.
    #[arg(short = 'a', long = "adapter", required = true)]
    adapters: Vec<String>,

    /// Minimum partial identity for end adapter matches (0.0-1.0).
    /// Uses identity over aligned region only - more lenient for partial adapters.
    #[arg(long, default_value_t = 0.75)]
    end_threshold: f64,

    /// Minimum full identity for middle adapter matches (0.0-1.0).
    /// Uses identity vs full adapter length - stricter than end matching.
    #[arg(long, default_value_t = 0.90)]
    middle_threshold: f64,

    /// Size of region to search at read ends (bp).
    #[arg(short = 'e', long, default_value_t = 150)]
    end_size: usize,

    /// Minimum alignment length to trigger trimming (bp).
    #[arg(long, default_value_t = 4)]
    min_trim_size: usize,

    /// Minimum fragment size after splitting (bp).
    #[arg(short = 'm', long, default_value_t = 1000)]
    min_length: usize,

    /// Extra bases to trim at read ends.
    #[arg(long, default_value_t = 2)]
    extra_end_trim: usize,

    /// Extra trim on expected side of middle adapter.
    #[arg(long, default_value_t = 10)]
    extra_middle_good: usize,

    /// Extra trim on unexpected side of middle adapter.
    #[arg(long, default_value_t = 100)]
    extra_middle_bad: usize,

    /// Split reads at internal adapters (chimera detection).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    split: bool,

    /// Number of reads per processing chunk (controls memory usage).
    #[arg(long, default_value_t = 1000)]
    chunk_size: usize,

    /// Number of threads for parallel processing.
    #[arg(short = 'j', long)]
    threads: Option<usize>,

    /// Gzip compression level (0-9). Auto-detected from output extension if not specified.
    #[arg(short = 'c', long)]
    compression_level: Option<u32>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Configure thread pool
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    // Prepare adapters
    let adapters = prepare_adapters(&args.adapters);
    eprintln!(
        "Loaded {} adapter(s) (searching both forward and reverse complement)",
        adapters.len()
    );

    // Build processing config
    let config = ProcessingConfig {
        end_size: args.end_size,
        end_threshold: args.end_threshold,
        middle_threshold: args.middle_threshold,
        min_trim_size: args.min_trim_size,
        min_read_length: args.min_length,
        extra_end_trim: args.extra_end_trim,
        extra_middle_good: args.extra_middle_good,
        extra_middle_bad: args.extra_middle_bad,
        split_enabled: args.split,
    };

    // Create reader
    let mut reader = if let Some(ref input_path) = args.input {
        parse_fastx_file(input_path)?
    } else {
        parse_fastx_stdin()?
    };

    // Create writer
    let writer: Box<dyn Write + Send> = if let Some(ref output_path) = args.output {
        let output_file = File::create(output_path)?;

        // Auto-detect compression from extension
        let should_compress = match args.compression_level {
            Some(0) => false,
            Some(_) => true,
            None => output_path.ends_with(".gz"),
        };

        if should_compress {
            let level = args.compression_level.unwrap_or(1);
            let mut builder = ParCompressBuilder::<Gzip>::new()
                .compression_level(GzpCompression::new(level as u32));

            if let Some(threads) = args.threads {
                builder = builder.num_threads(threads)?;
            }

            let encoder = builder.from_writer(output_file);
            Box::new(BufWriter::new(encoder))
        } else {
            Box::new(BufWriter::new(output_file))
        }
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };

    // Process reads in chunks
    let mut writer = writer;
    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(args.chunk_size);
    let mut total_reads = 0usize;
    let mut total_output = 0usize;
    let mut total_trimmed = 0usize;
    let mut total_split = 0usize;
    let mut total_discarded = 0usize;

    while let Some(record) = reader.next() {
        let rec = record?;

        chunk.push(FastqRecord {
            id: rec.id().to_vec(),
            seq: rec.seq().to_vec(),
            qual: rec.qual().map(|q| q.to_vec()).unwrap_or_default(),
        });

        if chunk.len() >= args.chunk_size {
            let stats = process_and_write_chunk(&mut chunk, &mut writer, &adapters, &config)?;
            total_reads += stats.0;
            total_output += stats.1;
            total_trimmed += stats.2;
            total_split += stats.3;
            total_discarded += stats.4;
            chunk.clear();
        }
    }

    // Process remaining reads
    if !chunk.is_empty() {
        let stats = process_and_write_chunk(&mut chunk, &mut writer, &adapters, &config)?;
        total_reads += stats.0;
        total_output += stats.1;
        total_trimmed += stats.2;
        total_split += stats.3;
        total_discarded += stats.4;
    }

    // Flush writer
    writer.flush()?;

    // Print summary
    eprintln!("\nProcessing complete:");
    eprintln!("  Input reads:    {}", total_reads);
    eprintln!("  Output reads:   {}", total_output);
    eprintln!("  Trimmed:        {}", total_trimmed);
    eprintln!("  Split:          {}", total_split);
    eprintln!("  Discarded:      {}", total_discarded);

    Ok(())
}

/// Process a chunk of reads in parallel and write results
fn process_and_write_chunk(
    chunk: &mut Vec<FastqRecord>,
    writer: &mut Box<dyn Write + Send>,
    adapters: &[Adapter],
    config: &ProcessingConfig,
) -> Result<(usize, usize, usize, usize, usize), Box<dyn std::error::Error>> {
    let input_count = chunk.len();

    // Process in parallel
    let results: Vec<ProcessedRead> = chunk
        .par_iter()
        .map(|record| process_read(record.clone(), adapters, config))
        .collect();

    // Collect stats and write sequentially
    let mut output_count = 0usize;
    let mut trimmed_count = 0usize;
    let mut split_count = 0usize;
    let mut discarded_count = 0usize;

    for processed in results {
        if processed.stats.was_discarded {
            discarded_count += 1;
        } else if processed.stats.was_split {
            split_count += 1;
        } else if processed.stats.bases_trimmed > 0 {
            trimmed_count += 1;
        }

        for record in processed.records {
            output_count += 1;
            write_fastq_record(writer, &record)?;
        }
    }

    Ok((input_count, output_count, trimmed_count, split_count, discarded_count))
}

/// Write a FASTQ record to the output
fn write_fastq_record(
    writer: &mut Box<dyn Write + Send>,
    record: &FastqRecord,
) -> Result<(), Box<dyn std::error::Error>> {
    writer.write_all(b"@")?;
    writer.write_all(&record.id)?;
    writer.write_all(b"\n")?;
    writer.write_all(&record.seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&record.qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}
