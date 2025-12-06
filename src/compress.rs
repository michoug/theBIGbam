//! Adaptive Run-Length Encoding (RLE) compression for genomic feature data.
//!
//! Compresses signals by grouping consecutive positions with similar values into runs.
//! For curves: uses self-referential threshold (compress_ratio * run_value).
//! For bars: uses coverage-relative threshold (compress_ratio * coverage[i]).
//! Bar values below threshold are filtered out (only significant deviations saved).

use crate::types::PlotType;

// ============================================================================
// RUN STRUCTURE
// ============================================================================

/// A run of consecutive positions with similar values.
#[derive(Debug, Clone)]
pub struct Run {
    /// First position in the run (1-indexed for database compatibility)
    pub start_pos: i32,
    /// Last position in the run (1-indexed, inclusive)
    pub end_pos: i32,
    /// The value for this run
    pub value: f32,
}

// ============================================================================
// MAIN COMPRESSION FUNCTION
// ============================================================================

/// Compress signal using adaptive run-length encoding.
///
/// # Algorithm
/// 1. Start first run with position 0 and its value
/// 2. For each subsequent position:
///    - If |value - current_run_value| <= compress_ratio * |current_run_value|, extend run
///    - Otherwise, close current run and start new one
/// 3. Return runs as (start_pos, end_pos, value) tuples
///
/// # Parameters
/// - `values`: The signal to compress (e.g., coverage at each position)
/// - `plot_type`: Curve (self-referential) or Bars (needs coverage context)
/// - `compress_ratio`: Relative tolerance (e.g., 0.1 for 10% change threshold)
///
/// # Returns
/// Vector of Run structs with 1-indexed positions
pub fn compress_signal(
    values: &[f64],
    plot_type: PlotType,
    compress_ratio: f64,
) -> Vec<Run> {
    compress_signal_with_reference(values, None, plot_type, compress_ratio)
}

/// Compress signal using adaptive RLE with optional coverage reference.
///
/// For bar plots (insertions, deletions, mismatches, read starts/ends),
/// compression should be relative to coverage depth, not the feature's own values.
///
/// # Parameters
/// - `values`: The signal to compress
/// - `reference`: Optional reference signal (e.g., coverage for bar plots)
/// - `plot_type`: Curve or Bars
/// - `compress_ratio`: Relative tolerance (e.g., 0.1 for 10% change threshold)
///
/// # Algorithm for Bars with Reference
/// For bar features, the threshold is `compress_ratio * reference[i]`.
/// This means a position with 5 mismatches at 1000x coverage (0.5%) gets compressed,
/// while 5 mismatches at 50x coverage (10%) creates a new run.
///
/// Bar values below the threshold are filtered out (not saved) - only significant
/// deviations from zero are preserved.
///
/// # Example
/// ```
/// coverage = [1000, 1000, 1000, 50, 50, 1000, 1000]
/// mismatches = [5, 5, 5, 5, 5, 5, 5]
/// compress_ratio = 0.1
///
/// Position 0-2: mismatch=5, coverage=1000, 5 < 0.1*1000 → filtered out (below threshold)
/// Position 3-4: mismatch=5, coverage=50, 5 > 0.1*50 → saved as run (spike!)
/// Position 5-6: mismatch=5, coverage=1000, 5 < 0.1*1000 → filtered out
///
/// Runs: [(4,5,5)]  // Only the spike at low coverage is preserved
/// ```
pub fn compress_signal_with_reference(
    values: &[f64],
    reference: Option<&[f64]>,
    plot_type: PlotType,
    compress_ratio: f64,
) -> Vec<Run> {
    let n = values.len();
    if n == 0 {
        return Vec::new();
    }

    // For bars with reference, use coverage-relative compression
    let use_reference = matches!(plot_type, PlotType::Bars) && reference.is_some();

    let mut runs = Vec::new();
    
    // Start first run
    let mut run_start = 0;
    let mut run_value = values[0];
    let mut run_sum = values[0];
    let mut run_count = 1;

    for i in 1..n {
        let val = values[i];
        
        // Determine threshold based on plot type and reference availability
        let threshold = if use_reference {
            // For bars: threshold is compress_ratio * coverage[i]
            let ref_val = reference.unwrap()[i];
            if ref_val.abs() < 1e-9 {
                compress_ratio  // If no coverage, use absolute threshold
            } else {
                compress_ratio * ref_val.abs()
            }
        } else {
            // For curves: threshold is compress_ratio * current_run_value (self-referential)
            if run_value.abs() < 1e-9 {
                compress_ratio  // Absolute threshold for near-zero values
            } else {
                compress_ratio * run_value.abs()  // Relative threshold
            }
        };

        if (val - run_value).abs() <= threshold {
            // Extend current run
            run_sum += val;
            run_count += 1;
            run_value = run_sum / run_count as f64;  // Update to running average
        } else {
            // Close current run and save it (if significant for bars)
            let is_significant = if use_reference {
                // For bars: only save if run_value is above threshold
                // (i.e., the value is significant relative to coverage)
                run_value.abs() > threshold
            } else {
                // For curves: always save
                true
            };

            if is_significant {
                runs.push(Run {
                    start_pos: (run_start + 1) as i32,  // Convert to 1-indexed
                    end_pos: (run_start + run_count) as i32,
                    value: run_value as f32,
                });
            }

            // Start new run
            run_start = i;
            run_value = val;
            run_sum = val;
            run_count = 1;
        }
    }

    // Don't forget to save the last run (if significant for bars)
    let final_threshold = if use_reference {
        let ref_val = reference.unwrap()[run_start];
        if ref_val.abs() < 1e-9 {
            compress_ratio
        } else {
            compress_ratio * ref_val.abs()
        }
    } else {
        if run_value.abs() < 1e-9 {
            compress_ratio
        } else {
            compress_ratio * run_value.abs()
        }
    };

    let is_significant = if use_reference {
        run_value.abs() > final_threshold
    } else {
        true
    };

    if is_significant {
        runs.push(Run {
            start_pos: (run_start + 1) as i32,
            end_pos: (run_start + run_count) as i32,
            value: run_value as f32,
        });
    }

    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adaptive_rle_simple() {
        // Test case from the documentation example
        let values = vec![100.0, 102.0, 98.0, 200.0, 205.0, 95.0, 100.0];
        let compress_ratio = 0.1; // 10% tolerance
        
        let runs = compress_signal(&values, PlotType::Curve, compress_ratio);
        
        // Expected runs (approximately):
        // Run 1: positions 1-3, value ~100
        // Run 2: positions 4-5, value ~202.5
        // Run 3: positions 6-7, value ~97.5
        assert_eq!(runs.len(), 3);
        
        assert_eq!(runs[0].start_pos, 1);
        assert_eq!(runs[0].end_pos, 3);
        assert!((runs[0].value - 100.0).abs() < 5.0); // ~100
        
        assert_eq!(runs[1].start_pos, 4);
        assert_eq!(runs[1].end_pos, 5);
        assert!((runs[1].value - 202.5).abs() < 5.0); // ~202.5
        
        assert_eq!(runs[2].start_pos, 6);
        assert_eq!(runs[2].end_pos, 7);
        assert!((runs[2].value - 97.5).abs() < 5.0); // ~97.5
    }

    #[test]
    fn test_adaptive_rle_constant() {
        // All values the same - should produce one run
        let values = vec![50.0; 1000];
        let compress_ratio = 0.1;
        
        let runs = compress_signal(&values, PlotType::Curve, compress_ratio);
        
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].start_pos, 1);
        assert_eq!(runs[0].end_pos, 1000);
        assert_eq!(runs[0].value, 50.0);
    }

    #[test]
    fn test_adaptive_rle_spike() {
        // Constant with a spike - spike should be its own run
        let values = vec![100.0, 100.0, 500.0, 100.0, 100.0];
        let compress_ratio = 0.1;
        
        let runs = compress_signal(&values, PlotType::Curve, compress_ratio);
        
        // Should have 3 runs: [1-2], [3], [4-5]
        assert_eq!(runs.len(), 3);
        
        assert_eq!(runs[0].start_pos, 1);
        assert_eq!(runs[0].end_pos, 2);
        
        assert_eq!(runs[1].start_pos, 3);
        assert_eq!(runs[1].end_pos, 3); // Singleton run
        assert_eq!(runs[1].value, 500.0);
        
        assert_eq!(runs[2].start_pos, 4);
        assert_eq!(runs[2].end_pos, 5);
    }

    #[test]
    fn test_adaptive_rle_empty() {
        let values: Vec<f64> = vec![];
        let compress_ratio = 0.1;
        
        let runs = compress_signal(&values, PlotType::Curve, compress_ratio);
        assert_eq!(runs.len(), 0);
    }

    #[test]
    fn test_context_aware_bars() {
        // Test bar compression with coverage reference and filtering
        // Scenario: constant mismatches (5) but varying coverage
        let coverage = vec![1000.0, 1000.0, 1000.0, 50.0, 50.0, 1000.0, 1000.0];
        let mismatches = vec![5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0];
        let compress_ratio = 0.1; // 10% threshold
        
        // With coverage reference:
        // Positions 0-2: 5 < 0.1*1000 → filtered out (insignificant at high coverage)
        // Position 3-4: 5 > 0.1*50 → saved as run (significant spike at low coverage!)
        // Position 5-6: 5 < 0.1*1000 → filtered out (back to insignificant)
        
        let runs = compress_signal_with_reference(&mismatches, Some(&coverage), PlotType::Bars, compress_ratio);
        
        // Should have only 1 run: the low-coverage spike (others filtered out)
        assert_eq!(runs.len(), 1);
        
        assert_eq!(runs[0].start_pos, 4);
        assert_eq!(runs[0].end_pos, 5); // Low coverage spike preserved
        assert_eq!(runs[0].value, 5.0);
    }

    #[test]
    fn test_bars_without_reference_fallback() {
        // Test that bars without reference still work (self-referential)
        let values = vec![5.0, 5.0, 50.0, 50.0, 5.0];
        let compress_ratio = 0.1;
        
        let runs = compress_signal(&values, PlotType::Bars, compress_ratio);
        
        // Should create runs based on value changes (10x jump)
        assert!(runs.len() >= 2); // At least spike separated
    }
}
