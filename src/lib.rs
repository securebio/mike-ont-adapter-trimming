//! ONT adapter trimming library
//!
//! Core algorithms for finding and trimming adapters from Oxford Nanopore reads.

use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alphabets::dna;

// ============================================================================
// Configuration
// ============================================================================

/// Porechop default scoring parameters
pub const MATCH_SCORE: i32 = 3;
pub const MISMATCH_SCORE: i32 = -6;
pub const GAP_OPEN: i32 = -5;
pub const GAP_EXTEND: i32 = -2;

/// Configuration for read processing
#[derive(Clone, Debug)]
pub struct ProcessingConfig {
    /// Size of region to search at read ends (default: 150bp)
    pub end_size: usize,
    /// Minimum partial identity for end adapter matches (default: 0.75)
    pub end_threshold: f64,
    /// Minimum full identity for middle adapter matches (default: 0.90)
    pub middle_threshold: f64,
    /// Minimum alignment length to trigger trimming (default: 4bp)
    pub min_trim_size: usize,
    /// Minimum fragment size after splitting (default: 1000bp)
    pub min_read_length: usize,
    /// Extra bases to trim at read ends (default: 2bp)
    pub extra_end_trim: usize,
    /// Extra trim on expected side of middle adapter (default: 10bp)
    pub extra_middle_good: usize,
    /// Extra trim on unexpected side of middle adapter (default: 100bp)
    pub extra_middle_bad: usize,
    /// Whether to split reads at internal adapters
    pub split_enabled: bool,
}

impl Default for ProcessingConfig {
    fn default() -> Self {
        Self {
            end_size: 150,
            end_threshold: 0.75,
            middle_threshold: 0.90,
            min_trim_size: 4,
            min_read_length: 1000,
            extra_end_trim: 2,
            extra_middle_good: 10,
            extra_middle_bad: 100,
            split_enabled: true,
        }
    }
}

// ============================================================================
// Data Structures
// ============================================================================

/// Adapter sequence with precomputed reverse complement
#[derive(Clone, Debug)]
pub struct Adapter {
    /// Name/identifier for this adapter
    pub name: String,
    /// Forward sequence
    pub seq: Vec<u8>,
    /// Reverse complement (precomputed)
    pub rev_comp: Vec<u8>,
}

impl Adapter {
    /// Create a new adapter with automatic reverse complement generation
    pub fn new(name: String, seq: Vec<u8>) -> Self {
        let rev_comp = dna::revcomp(&seq);
        Self { name, seq, rev_comp }
    }

    /// Create adapter from a sequence string
    pub fn from_str(name: &str, seq: &str) -> Self {
        Self::new(name.to_string(), seq.as_bytes().to_vec())
    }
}

/// Prepare adapters from CLI input strings
pub fn prepare_adapters(adapter_seqs: &[String]) -> Vec<Adapter> {
    adapter_seqs
        .iter()
        .enumerate()
        .map(|(i, seq)| Adapter::from_str(&format!("adapter_{}", i + 1), seq))
        .collect()
}

/// FASTQ record for processing
#[derive(Clone, Debug)]
pub struct FastqRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

/// Result of aligning an adapter to a sequence region
#[derive(Clone, Debug)]
pub struct AlignmentMatch {
    /// Index of the adapter that matched
    pub adapter_idx: usize,
    /// Start position in the read
    pub start: usize,
    /// End position in the read (exclusive)
    pub end: usize,
    /// Identity score (0.0 - 1.0)
    pub identity: f64,
    /// Whether the reverse complement matched
    pub is_reverse: bool,
    /// Aligned portion of adapter (start, end) for partial identity calculation
    pub adapter_aligned: (usize, usize),
}

/// Which end of the read we're searching
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EndRegion {
    Start,
    End,
}

/// Decision for how to process a read
#[derive(Clone, Debug)]
pub enum ProcessingDecision {
    /// No adapters found - keep read as-is
    NoTrim,
    /// Trim from start and/or end only
    TrimEnds {
        trim_start: usize,
        trim_end: usize,
    },
    /// Split read into multiple parts (adapter in middle)
    Split {
        /// (start, end) positions of each adapter to remove
        split_points: Vec<(usize, usize)>,
        trim_start: usize,
        trim_end: usize,
    },
    /// Discard entire read
    Discard {
        reason: String,
    },
}

/// Statistics from processing a read
#[derive(Clone, Debug, Default)]
pub struct ProcessingStats {
    pub adapters_found: usize,
    pub bases_trimmed: usize,
    pub was_split: bool,
    pub was_discarded: bool,
}

/// Result of processing a single read
#[derive(Clone, Debug)]
pub struct ProcessedRead {
    /// Output records (may be 0, 1, or multiple if split)
    pub records: Vec<FastqRecord>,
    /// Processing statistics
    pub stats: ProcessingStats,
}

// ============================================================================
// Alignment Functions
// ============================================================================

/// Result from a single alignment operation
#[derive(Clone, Debug)]
pub struct AlignmentResult {
    /// Alignment score
    pub score: i32,
    /// Start position in the read (x)
    pub read_start: usize,
    /// End position in the read (x, exclusive)
    pub read_end: usize,
    /// Start position in the adapter (y)
    pub adapter_start: usize,
    /// End position in the adapter (y, exclusive)
    pub adapter_end: usize,
}

/// Align an adapter to a sequence using custom fitting alignment
///
/// Uses custom alignment with free clipping at all four ends (matching Porechop ABI).
/// This allows partial adapters to be detected at read boundaries without penalty.
pub fn align_adapter(sequence: &[u8], adapter: &[u8]) -> AlignmentResult {
    // Create scoring with free clipping at all four ends (matching Porechop ABI)
    // xclip(0) = adapter can be partially aligned (free gaps at adapter ends)
    // yclip(0) = adapter can match anywhere in sequence (free gaps at sequence ends)
    let scoring = Scoring::from_scores(GAP_OPEN, GAP_EXTEND, MATCH_SCORE, MISMATCH_SCORE)
        .xclip(0)
        .yclip(0);

    let mut aligner = Aligner::with_scoring(scoring);

    // Custom alignment with adapter as x, sequence as y
    let alignment = aligner.custom(adapter, sequence);

    // Note: coordinates are swapped because adapter is x and sequence is y
    AlignmentResult {
        score: alignment.score,
        read_start: alignment.ystart,
        read_end: alignment.yend,
        adapter_start: alignment.xstart,
        adapter_end: alignment.xend,
    }
}

/// Calculate partial identity (identity over aligned region only)
/// Used for end adapter detection - more lenient for partial adapters
pub fn calculate_partial_identity(result: &AlignmentResult, _adapter_len: usize) -> f64 {
    let aligned_len = result.adapter_end - result.adapter_start;
    if aligned_len == 0 {
        return 0.0;
    }

    // Calculate matches from score
    // score = matches * MATCH_SCORE + mismatches * MISMATCH_SCORE + gaps * GAP_PENALTY
    // For a rough estimate, assume all aligned bases that aren't penalized are matches
    let max_score = aligned_len as i32 * MATCH_SCORE;
    let identity = if max_score > 0 {
        (result.score as f64 / max_score as f64).clamp(0.0, 1.0)
    } else {
        0.0
    };

    identity
}

/// Calculate full identity (identity vs full adapter length)
/// Used for middle adapter detection - stricter requirement
pub fn calculate_full_identity(result: &AlignmentResult, adapter_len: usize) -> f64 {
    if adapter_len == 0 {
        return 0.0;
    }

    let max_score = adapter_len as i32 * MATCH_SCORE;
    let identity = if max_score > 0 {
        (result.score as f64 / max_score as f64).clamp(0.0, 1.0)
    } else {
        0.0
    };

    identity
}

/// Find adapter match at read end (start or end region)
///
/// Uses partial identity scoring - more lenient for partial adapters at read ends.
pub fn find_adapter_at_end(
    sequence: &[u8],
    adapters: &[Adapter],
    region: EndRegion,
    config: &ProcessingConfig,
) -> Option<AlignmentMatch> {
    let seq_len = sequence.len();
    if seq_len == 0 {
        return None;
    }

    // Extract the region to search
    let (region_start, region_end) = match region {
        EndRegion::Start => (0, config.end_size.min(seq_len)),
        EndRegion::End => (seq_len.saturating_sub(config.end_size), seq_len),
    };

    let region_seq = &sequence[region_start..region_end];
    let region_len = region_end - region_start;

    let mut best_match: Option<AlignmentMatch> = None;
    let mut best_identity = config.end_threshold;

    for (adapter_idx, adapter) in adapters.iter().enumerate() {
        // Try forward orientation
        let result = align_adapter(region_seq, &adapter.seq);
        let adapter_aligned_len = result.adapter_end - result.adapter_start;
        let identity = calculate_partial_identity(&result, adapter.seq.len());

        // For end adapters, check:
        // 1. Identity above threshold
        // 2. Minimum alignment length
        // 3. Porechop-style boundary check: reject matches that span the entire search region
        //    (suspiciously long matches are likely false positives)
        //    - For Start region: match must NOT reach the end of search region
        //    - For End region: match must NOT start at beginning of search region
        let is_valid_match = match region {
            EndRegion::Start => result.read_end != region_len,
            EndRegion::End => result.read_start != 0,
        };

        if identity > best_identity
            && adapter_aligned_len >= config.min_trim_size
            && is_valid_match
        {
            best_identity = identity;
            best_match = Some(AlignmentMatch {
                adapter_idx,
                start: region_start + result.read_start,
                end: region_start + result.read_end,
                identity,
                is_reverse: false,
                adapter_aligned: (result.adapter_start, result.adapter_end),
            });
        }

        // Try reverse complement
        let result_rc = align_adapter(region_seq, &adapter.rev_comp);
        let adapter_aligned_len_rc = result_rc.adapter_end - result_rc.adapter_start;
        let identity_rc = calculate_partial_identity(&result_rc, adapter.rev_comp.len());

        let is_valid_match_rc = match region {
            EndRegion::Start => result_rc.read_end != region_len,
            EndRegion::End => result_rc.read_start != 0,
        };

        if identity_rc > best_identity
            && adapter_aligned_len_rc >= config.min_trim_size
            && is_valid_match_rc
        {
            best_identity = identity_rc;
            best_match = Some(AlignmentMatch {
                adapter_idx,
                start: region_start + result_rc.read_start,
                end: region_start + result_rc.read_end,
                identity: identity_rc,
                is_reverse: true,
                adapter_aligned: (result_rc.adapter_start, result_rc.adapter_end),
            });
        }
    }

    best_match
}

/// Find all internal adapter matches (for chimera detection/splitting)
///
/// Uses full identity scoring - stricter requirement for middle adapters.
/// Uses iterative masking to find multiple adapters in same read.
pub fn find_adapters_internal(
    sequence: &[u8],
    adapters: &[Adapter],
    config: &ProcessingConfig,
) -> Vec<AlignmentMatch> {
    let seq_len = sequence.len();
    if seq_len <= 2 * config.end_size {
        return Vec::new(); // Read too short to have middle region
    }

    let middle_start = config.end_size;
    let middle_end = seq_len - config.end_size;

    // Work with a mutable copy for masking
    let mut masked_seq = sequence.to_vec();
    let mut matches = Vec::new();

    loop {
        let mut best_match: Option<AlignmentMatch> = None;
        let mut best_identity = config.middle_threshold;

        // Search middle region
        let middle_region = &masked_seq[middle_start..middle_end];

        for (adapter_idx, adapter) in adapters.iter().enumerate() {
            // Try forward orientation
            let result = align_adapter(middle_region, &adapter.seq);
            let identity = calculate_full_identity(&result, adapter.seq.len());

            if identity > best_identity {
                best_identity = identity;
                best_match = Some(AlignmentMatch {
                    adapter_idx,
                    start: middle_start + result.read_start,
                    end: middle_start + result.read_end,
                    identity,
                    is_reverse: false,
                    adapter_aligned: (result.adapter_start, result.adapter_end),
                });
            }

            // Try reverse complement
            let result_rc = align_adapter(middle_region, &adapter.rev_comp);
            let identity_rc = calculate_full_identity(&result_rc, adapter.rev_comp.len());

            if identity_rc > best_identity {
                best_identity = identity_rc;
                best_match = Some(AlignmentMatch {
                    adapter_idx,
                    start: middle_start + result_rc.read_start,
                    end: middle_start + result_rc.read_end,
                    identity: identity_rc,
                    is_reverse: true,
                    adapter_aligned: (result_rc.adapter_start, result_rc.adapter_end),
                });
            }
        }

        match best_match {
            Some(m) => {
                // Mask the found adapter region with Ns
                for i in m.start..m.end {
                    masked_seq[i] = b'N';
                }
                matches.push(m);
            }
            None => break, // No more matches found
        }
    }

    // Sort by position
    matches.sort_by_key(|m| m.start);
    matches
}

// ============================================================================
// Processing Functions
// ============================================================================

/// Scan a read for all adapter matches
pub fn scan_read_for_adapters(
    sequence: &[u8],
    adapters: &[Adapter],
    config: &ProcessingConfig,
) -> Vec<AlignmentMatch> {
    let mut all_matches = Vec::new();

    // Check start region
    if let Some(m) = find_adapter_at_end(sequence, adapters, EndRegion::Start, config) {
        all_matches.push(m);
    }

    // Check end region
    if let Some(m) = find_adapter_at_end(sequence, adapters, EndRegion::End, config) {
        all_matches.push(m);
    }

    // Check middle region for splitting
    if config.split_enabled {
        let middle_matches = find_adapters_internal(sequence, adapters, config);
        all_matches.extend(middle_matches);
    }

    all_matches
}

/// Decide how to process a read based on adapter matches
pub fn decide_processing(
    matches: &[AlignmentMatch],
    read_length: usize,
    config: &ProcessingConfig,
) -> ProcessingDecision {
    if matches.is_empty() {
        return ProcessingDecision::NoTrim;
    }

    let mut trim_start = 0usize;
    let mut trim_end = 0usize;
    let mut split_points: Vec<(usize, usize)> = Vec::new();

    for m in matches {
        // Categorize by position
        let is_at_start = m.end <= config.end_size;
        let is_at_end = m.start >= read_length.saturating_sub(config.end_size);
        let is_middle = !is_at_start && !is_at_end;

        if is_at_start {
            // Trim from start: adapter end + extra trim
            let new_trim = m.end + config.extra_end_trim;
            trim_start = trim_start.max(new_trim);
        } else if is_at_end {
            // Trim from end: everything from adapter start - extra trim
            let trim_from_end = (read_length - m.start) + config.extra_end_trim;
            trim_end = trim_end.max(trim_from_end);
        } else if is_middle {
            // Add to split points with extra trimming
            // For middle adapters, we trim extra on both sides
            let split_start = m.start.saturating_sub(config.extra_middle_bad);
            let split_end = (m.end + config.extra_middle_good).min(read_length);
            split_points.push((split_start, split_end));
        }
    }

    // Check if anything remains after trimming
    let remaining = read_length.saturating_sub(trim_start).saturating_sub(trim_end);
    if remaining < config.min_read_length {
        return ProcessingDecision::Discard {
            reason: "Read too short after trimming".to_string(),
        };
    }

    if !split_points.is_empty() {
        ProcessingDecision::Split {
            split_points,
            trim_start,
            trim_end,
        }
    } else if trim_start > 0 || trim_end > 0 {
        ProcessingDecision::TrimEnds {
            trim_start,
            trim_end,
        }
    } else {
        ProcessingDecision::NoTrim
    }
}

/// Process a single read: find adapters, trim/split, generate output records
pub fn process_read(
    record: FastqRecord,
    adapters: &[Adapter],
    config: &ProcessingConfig,
) -> ProcessedRead {
    let seq_len = record.seq.len();
    let matches = scan_read_for_adapters(&record.seq, adapters, config);
    let decision = decide_processing(&matches, seq_len, config);

    let mut stats = ProcessingStats {
        adapters_found: matches.len(),
        ..Default::default()
    };

    let records = match decision {
        ProcessingDecision::NoTrim => {
            vec![record]
        }
        ProcessingDecision::TrimEnds { trim_start, trim_end } => {
            let new_end = seq_len.saturating_sub(trim_end);
            if trim_start >= new_end {
                stats.was_discarded = true;
                Vec::new()
            } else {
                stats.bases_trimmed = trim_start + trim_end;
                vec![FastqRecord {
                    id: record.id,
                    seq: record.seq[trim_start..new_end].to_vec(),
                    qual: record.qual[trim_start..new_end].to_vec(),
                }]
            }
        }
        ProcessingDecision::Split { split_points, trim_start, trim_end } => {
            stats.was_split = true;

            // Build list of segments to keep
            let mut segments: Vec<(usize, usize)> = Vec::new();
            let new_end = seq_len.saturating_sub(trim_end);

            // Start with the region after initial trim
            let mut current_start = trim_start;

            for (split_start, split_end) in &split_points {
                if current_start < *split_start {
                    segments.push((current_start, *split_start));
                }
                current_start = *split_end;
            }

            // Add final segment
            if current_start < new_end {
                segments.push((current_start, new_end));
            }

            // Filter by minimum length and generate records
            let mut output_records = Vec::new();
            let base_id = String::from_utf8_lossy(&record.id);

            for (i, (start, end)) in segments.iter().enumerate() {
                let len = end - start;
                if len >= config.min_read_length {
                    let new_id = if segments.len() > 1 {
                        format!("{}_part{}", base_id, i + 1)
                    } else {
                        base_id.to_string()
                    };

                    output_records.push(FastqRecord {
                        id: new_id.into_bytes(),
                        seq: record.seq[*start..*end].to_vec(),
                        qual: record.qual[*start..*end].to_vec(),
                    });
                    stats.bases_trimmed += seq_len - len;
                }
            }

            if output_records.is_empty() {
                stats.was_discarded = true;
            }

            output_records
        }
        ProcessingDecision::Discard { reason: _ } => {
            stats.was_discarded = true;
            Vec::new()
        }
    };

    ProcessedRead { records, stats }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adapter_reverse_complement() {
        let adapter = Adapter::from_str("test", "ACGT");
        assert_eq!(adapter.seq, b"ACGT");
        assert_eq!(adapter.rev_comp, b"ACGT"); // ACGT is its own reverse complement

        let adapter2 = Adapter::from_str("test2", "AACG");
        assert_eq!(adapter2.rev_comp, b"CGTT");
    }

    #[test]
    fn test_prepare_adapters() {
        let seqs = vec!["ACGT".to_string(), "AAAA".to_string()];
        let adapters = prepare_adapters(&seqs);
        assert_eq!(adapters.len(), 2);
        assert_eq!(adapters[0].name, "adapter_1");
        assert_eq!(adapters[1].name, "adapter_2");
    }

    #[test]
    fn test_no_trim_when_no_adapter() {
        let config = ProcessingConfig::default();
        let adapters = vec![Adapter::from_str("test", "NNNNNNNNNN")]; // Won't match

        let record = FastqRecord {
            id: b"read1".to_vec(),
            seq: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            qual: vec![b'I'; 200],
        };

        let result = process_read(record.clone(), &adapters, &config);
        assert_eq!(result.records.len(), 1);
        assert_eq!(result.records[0].seq, record.seq);
    }
}
