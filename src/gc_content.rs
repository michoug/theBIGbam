//! GC content computation for genomic sequences.
//!
//! Computes GC content (percentage of G and C bases) using a sliding window approach
//! and applies RLE compression for efficient storage.

/// A run of consecutive positions with similar GC percentage.
#[derive(Debug, Clone)]
pub struct GCContentRun {
    /// First position in the run (1-indexed for database compatibility)
    pub start_pos: i32,
    /// Last position in the run (1-indexed, inclusive)
    pub end_pos: i32,
    /// GC percentage for this run (0-100)
    pub gc_percentage: u8,
}

/// Statistics for GC content across a contig.
#[derive(Debug, Clone)]
pub struct GCStats {
    /// Mean GC percentage (0-100)
    pub average: f32,
    /// Standard deviation of windowed GC percentages
    pub sd: f32,
    /// Median GC percentage (0-100)
    pub median: f32,
}

/// Compute GC content using a sliding window approach with RLE compression.
///
/// # Parameters
/// - `sequence`: The DNA sequence as bytes (A, T, G, C, N)
/// - `window_size`: Size of the sliding window in base pairs (typically 100)
/// - `contig_variation_percentage`: RLE compression tolerance for contig-level features (default 0.1%)
///
/// # Algorithm
/// 1. For each position, compute GC% in a centered window
/// 2. N bases are excluded from both numerator and denominator
/// 3. Apply RLE compression: merge consecutive positions with similar GC%
///
/// # Returns
/// Tuple of (runs, stats) where:
/// - runs: Vector of GC content runs with (start_pos, end_pos, gc_percentage)
/// - stats: GCStats with average, sd, and median GC percentages
pub fn compute_gc_content(sequence: &[u8], window_size: usize, contig_variation_percentage: f64) -> (Vec<GCContentRun>, GCStats) {
    let n = sequence.len();
    if n == 0 {
        return (Vec::new(), GCStats { average: 0.0, sd: 0.0, median: 0.0 });
    }

    let half_window = window_size / 2;
    let mut gc_values: Vec<u8> = Vec::with_capacity(n);

    // Compute GC content for each position using a centered sliding window
    for i in 0..n {
        let start = i.saturating_sub(half_window);
        let end = (i + half_window + 1).min(n);

        let (gc_count, valid_count) = count_gc_in_window(&sequence[start..end]);

        let gc_pct = if valid_count > 0 {
            ((gc_count as f64 / valid_count as f64) * 100.0).round() as u8
        } else {
            50 // Default to 50% if window contains only Ns
        };

        gc_values.push(gc_pct);
    }

    // Compute statistics from raw values before compression
    let stats = compute_gc_stats(&gc_values);

    // Apply RLE compression with user-configurable tolerance
    let runs = compress_gc_values(&gc_values, contig_variation_percentage);

    (runs, stats)
}

/// Compute GC statistics (average, sd, median) from raw GC values.
fn compute_gc_stats(gc_values: &[u8]) -> GCStats {
    if gc_values.is_empty() {
        return GCStats { average: 0.0, sd: 0.0, median: 0.0 };
    }

    let n = gc_values.len() as f64;

    // Compute average
    let sum: f64 = gc_values.iter().map(|&v| v as f64).sum();
    let average = sum / n;

    // Compute standard deviation
    let variance: f64 = gc_values.iter()
        .map(|&v| {
            let diff = v as f64 - average;
            diff * diff
        })
        .sum::<f64>() / n;
    let sd = variance.sqrt();

    // Compute median
    let mut sorted = gc_values.to_vec();
    sorted.sort_unstable();
    let median = if sorted.len() % 2 == 0 {
        let mid = sorted.len() / 2;
        (sorted[mid - 1] as f64 + sorted[mid] as f64) / 2.0
    } else {
        sorted[sorted.len() / 2] as f64
    };

    GCStats {
        average: average as f32,
        sd: sd as f32,
        median: median as f32,
    }
}

/// Count G and C bases in a window, excluding N bases.
#[inline]
fn count_gc_in_window(window: &[u8]) -> (usize, usize) {
    let mut gc_count = 0;
    let mut valid_count = 0;

    for &base in window {
        match base {
            b'G' | b'g' | b'C' | b'c' => {
                gc_count += 1;
                valid_count += 1;
            }
            b'A' | b'a' | b'T' | b't' => {
                valid_count += 1;
            }
            // N and other characters are excluded
            _ => {}
        }
    }

    (gc_count, valid_count)
}

/// Compress GC values using adaptive RLE.
/// Uses the formula: |x[i] - x[i-1]| <= ratio × min(x[i], x[i-1])
/// Input is percentage (e.g., 10 for 10%), converted to fraction internally.
fn compress_gc_values(values: &[u8], percentage: f64) -> Vec<GCContentRun> {
    if values.is_empty() {
        return Vec::new();
    }

    let ratio = percentage * 0.01; // Convert percentage to fraction
    let mut runs = Vec::new();

    let mut run_start = 0;
    let mut run_value = values[0] as f64;
    let mut run_sum = values[0] as f64;
    let mut run_count = 1;

    for i in 1..values.len() {
        let val = values[i] as f64;

        // RLE formula: |x[i] - x[i-1]| <= ratio × min(x[i], x[i-1])
        let min_val = val.min(run_value);
        let threshold = ratio * min_val.max(1.0); // Use at least 1.0 to handle zero values

        if (val - run_value).abs() <= threshold {
            // Extend current run
            run_sum += val;
            run_count += 1;
            run_value = run_sum / run_count as f64;
        } else {
            // Close current run and save it
            runs.push(GCContentRun {
                start_pos: (run_start + 1) as i32,
                end_pos: (run_start + run_count) as i32,
                gc_percentage: run_value.round() as u8,
            });

            // Start new run
            run_start = i;
            run_value = val;
            run_sum = val;
            run_count = 1;
        }
    }

    // Save the last run
    runs.push(GCContentRun {
        start_pos: (run_start + 1) as i32,
        end_pos: (run_start + run_count) as i32,
        gc_percentage: run_value.round() as u8,
    });

    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content_simple() {
        // ATGC sequence - 50% GC
        let sequence = b"ATGCATGCATGC";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with ~50% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert!(run.gc_percentage >= 40 && run.gc_percentage <= 60);
        }
        // Stats should also be around 50%
        assert!(stats.average >= 40.0 && stats.average <= 60.0);
    }

    #[test]
    fn test_gc_content_high_gc() {
        // All GC
        let sequence = b"GGGGCCCCGGGGCCCC";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with 100% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert_eq!(run.gc_percentage, 100);
        }
        // Stats should be 100% with 0 sd
        assert_eq!(stats.average, 100.0);
        assert_eq!(stats.sd, 0.0);
        assert_eq!(stats.median, 100.0);
    }

    #[test]
    fn test_gc_content_low_gc() {
        // All AT
        let sequence = b"AAAATTTTAAAATTTT";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);

        // Should produce runs with 0% GC
        assert!(!runs.is_empty());
        for run in &runs {
            assert_eq!(run.gc_percentage, 0);
        }
        // Stats should be 0% with 0 sd
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
        assert_eq!(stats.median, 0.0);
    }

    #[test]
    fn test_gc_content_with_n() {
        // Sequence with N bases - N should be excluded
        let sequence = b"ATGCNNNNATGC";
        let (runs, _stats) = compute_gc_content(sequence, 4, 10.0);

        // Should still compute ~50% GC from valid bases
        assert!(!runs.is_empty());
    }

    #[test]
    fn test_gc_content_empty() {
        let sequence: &[u8] = b"";
        let (runs, stats) = compute_gc_content(sequence, 4, 10.0);
        assert!(runs.is_empty());
        assert_eq!(stats.average, 0.0);
        assert_eq!(stats.sd, 0.0);
        assert_eq!(stats.median, 0.0);
    }

    #[test]
    fn test_full_coverage() {
        // Verify that the runs cover the entire sequence
        let sequence = b"ATGCATGCATGCATGCATGC"; // 20 bases
        let (runs, _stats) = compute_gc_content(sequence, 4, 10.0);

        // First run should start at position 1
        assert_eq!(runs.first().unwrap().start_pos, 1);
        // Last run should end at position 20
        assert_eq!(runs.last().unwrap().end_pos, 20);

        // Verify continuous coverage (no gaps)
        for i in 1..runs.len() {
            assert_eq!(runs[i].start_pos, runs[i - 1].end_pos + 1);
        }
    }

    #[test]
    fn test_compression() {
        // Long sequence with uniform GC should compress well
        let sequence: Vec<u8> = b"ATGC".repeat(1000);
        let (runs, stats) = compute_gc_content(&sequence, 100, 10.0);

        // Should compress but still cover full length
        assert!(!runs.is_empty());
        assert_eq!(runs.first().unwrap().start_pos, 1);
        assert_eq!(runs.last().unwrap().end_pos, 4000);
        // Stats should be around 50%
        assert!(stats.average >= 45.0 && stats.average <= 55.0);
    }

    #[test]
    fn test_gc_stats() {
        // Test stats computation with known values
        let gc_values: Vec<u8> = vec![40, 50, 60];
        let stats = compute_gc_stats(&gc_values);

        // Average should be 50
        assert_eq!(stats.average, 50.0);
        // Median should be 50
        assert_eq!(stats.median, 50.0);
        // SD should be sqrt((100+0+100)/3) = sqrt(200/3) ≈ 8.16
        assert!(stats.sd > 8.0 && stats.sd < 9.0);
    }
}
