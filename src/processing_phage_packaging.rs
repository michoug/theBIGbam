//! Phage packaging mechanism classification and terminal repeat detection.

use crate::db::RepeatsData;
use crate::types::TerminusArea;

/// Peak for phage termini detection (legacy, used internally).
#[derive(Clone, Debug)]
pub struct Peak {
    pub position: i32,
    pub value: f64,
}

/// Diagnostic information for a terminus position (for CSV export).
#[derive(Clone, Debug)]
pub struct TerminusDiagnostic {
    pub position: i32,
    pub terminus_type: String,  // "start" or "end"
    pub spc: u64,
    pub clippings: u64,
    pub observed_ratio: f64,
    pub primary_reads: u64,
    pub coverage_reduced: u64,
    pub tau: f64,
    pub kept: bool,
    pub discarded_because: String,
}

/// Peak area for phage termini detection.
/// Represents a merged region of nearby positions with significant signal.
#[derive(Clone, Debug)]
pub struct PeakArea {
    /// First position of the merged area
    pub start_pos: i32,
    /// Last position of the merged area
    pub end_pos: i32,
    /// Position with highest SPC (read starts + read ends)
    pub center_pos: i32,
    /// Total SPC (sum of read_starts + read_ends in area)
    pub total_spc: u32,
    /// Total clippings in area (for informational purposes)
    pub total_clips: u32,
    /// Coverage at center position
    pub coverage: u32,
    /// Tau value at center position
    pub tau: f64,
}

/// Configuration for phage termini detection.
#[derive(Clone, Copy)]
pub struct PhageTerminiConfig {
    /// Minimum aligned fraction to evaluate the phage termini
    pub min_aligned_fraction: i32,
    /// Minimum identity (%) for considering duplications as DTR/ITR
    pub min_identity_dtr: i32,
    /// Maximum distance (bp) from reference ends for duplication regions to be considered valid.
    /// Duplications must have one region starting within this distance from position 1
    /// and the other region ending within this distance from the contig end.
    pub max_distance_duplication: i32,
    /// Minimum event count to consider a signal significant (applies to both peaks and clippings)
    pub min_events: i32,
    /// Minimum frequence of events to consider a signal significant (applies to both peaks and clippings)
    pub min_frequency: i32,
    /// Maximum distance (bp) between peaks to merge them
    pub max_distance_peaks: i32,
    /// Significance threshold for clipping test (p-value threshold)
    /// Lower values are more stringent (require higher confidence that peak is real terminus)
    pub clipping_significance: f64,
}

impl Default for PhageTerminiConfig {
    fn default() -> Self {
        Self {
            min_aligned_fraction: 90,
            min_identity_dtr: 90,
            max_distance_duplication: 100,
            min_events: 10,
            min_frequency: 10,
            max_distance_peaks: 20,
            clipping_significance: 0.05, // 95% confidence interval
        }
    }
}

/// Calculate circular distance between two positions.
fn circular_distance(pos1: i32, pos2: i32, genome_length: usize) -> i32 {
    let gl = genome_length as i32;
    let dist_forward = if pos2 >= pos1 {
        pos2 - pos1
    } else {
        (gl - pos1) + pos2
    };
    let dist_backward = if pos1 >= pos2 {
        pos1 - pos2
    } else {
        (gl - pos2) + pos1
    };
    dist_forward.min(dist_backward)
}

/// DTR region info for distance calculations.
/// first_start/first_end: positions of the first (kept) duplicated region
/// second_start/second_end: positions of the second (zeroed) duplicated region
#[derive(Clone, Debug)]
pub struct DtrRegion {
    pub first_start: i32,
    pub first_end: i32,
    pub second_start: i32,
    pub second_end: i32,
    /// true = DTR (direct), false = ITR (inverted)
    pub is_direct: bool,
}

/// Calculate minimum distance between two positions, considering DTR equivalence.
///
/// If a position is near the second DTR region, it's treated as being near
/// the equivalent position in the first DTR region.
fn dtr_aware_distance(
    pos1: i32,
    pos2: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
    max_distance: i32,
) -> i32 {
    // Start with regular distance
    let regular_dist = if circular {
        circular_distance(pos1, pos2, genome_length)
    } else {
        (pos1 - pos2).abs()
    };

    // If no DTR regions or already close enough, return regular distance
    if dtr_regions.is_empty() || regular_dist <= max_distance {
        return regular_dist;
    }

    // Check if either position is near a second DTR region
    // If so, calculate distance using the equivalent first region position
    let mut min_dist = regular_dist;

    for dtr in dtr_regions {
        // Check if pos1 is near second region - calculate virtual position
        let pos1_virtual = get_virtual_position(pos1, dtr, max_distance);
        // Check if pos2 is near second region - calculate virtual position
        let pos2_virtual = get_virtual_position(pos2, dtr, max_distance);

        // Calculate distance using virtual positions
        let virtual_dist = if circular {
            circular_distance(pos1_virtual, pos2_virtual, genome_length)
        } else {
            (pos1_virtual - pos2_virtual).abs()
        };

        min_dist = min_dist.min(virtual_dist);
    }

    min_dist
}

/// Get the virtual position for a point near a terminal repeat region.
/// If the position is near the second region, return the equivalent
/// position in the first region. Otherwise, return the original position.
///
/// Mapping differs by repeat type:
/// - DTR (direct):   second_start + offset → first_start + offset
/// - ITR (inverted): second_start + offset → first_end - offset
fn get_virtual_position(pos: i32, dtr: &DtrRegion, max_distance: i32) -> i32 {
    // Check if position is in or near second region
    let in_second = pos >= dtr.second_start && pos <= dtr.second_end;
    let near_second_start = (pos - dtr.second_start).abs() <= max_distance;
    let near_second_end = (pos - dtr.second_end).abs() <= max_distance;

    if !in_second && !near_second_start && !near_second_end {
        // Position is not near second region, keep as is
        return pos;
    }

    // Calculate offset from second region start
    let offset = pos - dtr.second_start;

    if dtr.is_direct {
        // DTR: same direction mapping
        dtr.first_start + offset
    } else {
        // ITR: inverted mapping
        dtr.first_end - offset
    }
}

/// Statistical test for clipping artifacts using binomial proportion test.
///
/// Tests whether the observed ratio of read_starts/(read_starts+clippings)
/// is significantly LOWER than the expected unclipped_ratio (i.e., has
/// significantly more clippings than expected).
///
/// H0: ratio = unclipped_ratio (position is normal)
/// H1: ratio < unclipped_ratio (position has excess clippings → artifact)
///
/// Returns true if position should be KEPT (not a clipping artifact)
/// Returns false if position should be DISCARDED (likely clipping artifact)
fn is_significant_terminus(
    read_starts: u32,
    clippings: u32,
    unclipped_ratio: f64,
    clipping_significance: f64,
) -> bool {
    let n = read_starts + clippings;
    if n == 0 {
        return false;
    }
    if clippings == 0 {
        return true; // No clipping → keep
    }

    let p_obs = read_starts as f64 / n as f64;
    let p_exp = unclipped_ratio;

    // If observed is at or above expected, this position has same or fewer
    // clippings than average → definitely not a clipping artifact → keep
    if p_obs >= p_exp {
        return true;
    }

    // Observed is below expected - test if SIGNIFICANTLY more clippings
    // Standard error under null hypothesis
    let std_err = (p_exp * (1.0 - p_exp) / n as f64).sqrt();
    if std_err == 0.0 {
        return true; // Can't compute, keep by default
    }

    // Z-score: how many std errors below expected?
    let z = (p_exp - p_obs) / std_err;

    // Critical z-value for one-tailed test
    let z_critical = if clipping_significance >= 0.10 {
        1.282 // 90% confidence
    } else if clipping_significance >= 0.05 {
        1.645 // 95% confidence
    } else if clipping_significance >= 0.01 {
        2.326 // 99% confidence
    } else {
        2.576 // 99.5% confidence
    };

    // If z > z_critical, position has significantly MORE clippings → discard
    // Otherwise, not significantly worse → keep
    z <= z_critical
}

/// Statistical test for clipping artifacts with detailed reason for rejection.
/// Returns (passes_test, rejection_reason).
///
/// Tests if position has significantly MORE clippings than expected.
/// If yes → discard (clipping artifact). If no → keep.
fn is_significant_terminus_with_reason(
    read_starts: u32,
    clippings: u32,
    unclipped_ratio: f64,
    clipping_significance: f64,
) -> (bool, String) {
    let n = read_starts + clippings;
    if n == 0 {
        return (false, "no_data".to_string());
    }
    if clippings == 0 {
        return (true, String::new()); // No clipping → keep
    }

    let p_obs = read_starts as f64 / n as f64;
    let p_exp = unclipped_ratio;

    // If observed is at or above expected, this position has same or fewer
    // clippings than average → definitely not a clipping artifact → keep
    if p_obs >= p_exp {
        return (true, String::new());
    }

    // Observed is below expected - test if SIGNIFICANTLY more clippings
    // Standard error under null hypothesis
    let std_err = (p_exp * (1.0 - p_exp) / n as f64).sqrt();
    if std_err == 0.0 {
        return (true, String::new()); // Can't compute, keep by default
    }

    // Z-score: how many std errors below expected?
    let z = (p_exp - p_obs) / std_err;

    // Critical z-value for one-tailed test
    let z_critical = if clipping_significance >= 0.10 {
        1.282 // 90% confidence
    } else if clipping_significance >= 0.05 {
        1.645 // 95% confidence
    } else if clipping_significance >= 0.01 {
        2.326 // 99% confidence
    } else {
        2.576 // 99.5% confidence
    };

    // If z > z_critical, position has significantly MORE clippings → discard
    if z > z_critical {
        (false, format!("statistical_test: excess_clippings z={:.2} > z_critical={:.2}", z, z_critical))
    } else {
        (true, String::new()) // Not significantly worse → keep
    }
}

/// Sum clippings within max_distance of a center position.
/// This ensures that clipping events near a peak are properly associated with it.
fn sum_clippings_in_range(clippings: &[u64], center_idx: usize, max_distance: i32) -> u64 {
    let start = center_idx.saturating_sub(max_distance as usize);
    let end = (center_idx + max_distance as usize + 1).min(clippings.len());
    clippings[start..end].iter().sum()
}

/// Filter positions and merge nearby ones into peak areas.
///
/// This implements the new flow:
/// 1. Filter positions: read_starts >= 0.1*coverage_reduced AND read_starts >= min_events
/// 2. Filter positions: statistical clipping test (with clipping pre-filter)
/// 3. Merge nearby positions into peak-areas (DTR-aware)
/// 4. Filter areas: total_spc >= 0.1 * coverage_reduced_mean
///
/// Returns (start_areas, end_areas) for classification.
pub fn filter_and_merge_to_areas(
    reads_starts: &[u64],
    reads_ends: &[u64],
    left_clippings: &[u64],
    right_clippings: &[u64],
    coverage_reduced: &[u64],
    primary_reads: &[u64],
    unclipped_ratio: f64,
    _coverage_reduced_mean: f64,
    config: &PhageTerminiConfig,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (Vec<PeakArea>, Vec<PeakArea>) {
    let min_events = config.min_events as u64;
    let max_distance_peaks = config.max_distance_peaks;
    let clipping_significance = config.clipping_significance;

    // Step 1-2: Filter positions for starts
    let mut filtered_start_positions: Vec<(i32, u64, u64)> = Vec::new(); // (pos, spc, clips)
    for (idx, &spc) in reads_starts.iter().enumerate() {
        if spc == 0 {
            continue;
        }
        let pos = (idx + 1) as i32; // 1-indexed
        let cov = coverage_reduced.get(idx).copied().unwrap_or(0);
        // Sum clippings within max_distance_peaks of this position
        let clips = sum_clippings_in_range(left_clippings, idx, max_distance_peaks);
        let prim = primary_reads.get(idx).copied().unwrap_or(0);

        // Step 1: Local criteria
        let local_threshold = (cov as f64 * 0.1) as u64;
        if spc < local_threshold || spc < min_events {
            continue;
        }

        // Pre-filter clippings: only consider them significant if they meet threshold
        let effective_clips = if clips >= (prim as f64 * 0.1) as u64 && clips >= min_events {
            clips
        } else {
            0
        };

        // Step 2: Statistical test
        if !is_significant_terminus(spc as u32, effective_clips as u32, unclipped_ratio, clipping_significance) {
            continue;
        }

        filtered_start_positions.push((pos, spc, effective_clips));
    }

    // Step 1-2: Filter positions for ends
    let mut filtered_end_positions: Vec<(i32, u64, u64)> = Vec::new();
    for (idx, &spc) in reads_ends.iter().enumerate() {
        if spc == 0 {
            continue;
        }
        let pos = (idx + 1) as i32;
        let cov = coverage_reduced.get(idx).copied().unwrap_or(0);
        // Sum clippings within max_distance_peaks of this position
        let clips = sum_clippings_in_range(right_clippings, idx, max_distance_peaks);
        let prim = primary_reads.get(idx).copied().unwrap_or(0);

        let local_threshold = (cov as f64 * 0.1) as u64;
        if spc < local_threshold || spc < min_events {
            continue;
        }

        // Pre-filter clippings
        let effective_clips = if clips >= (prim as f64 * 0.1) as u64 && clips >= min_events {
            clips
        } else {
            0
        };

        if !is_significant_terminus(spc as u32, effective_clips as u32, unclipped_ratio, clipping_significance) {
            continue;
        }

        filtered_end_positions.push((pos, spc, effective_clips));
    }

    // Step 3: Merge nearby positions into areas (starts)
    let start_areas = merge_positions_to_areas(
        &filtered_start_positions,
        reads_starts,
        reads_ends,
        left_clippings,
        coverage_reduced,
        max_distance_peaks,
        genome_length,
        circular,
        dtr_regions,
    );

    // Step 3: Merge nearby positions into areas (ends)
    let end_areas = merge_positions_to_areas(
        &filtered_end_positions,
        reads_starts,
        reads_ends,
        right_clippings,
        coverage_reduced,
        max_distance_peaks,
        genome_length,
        circular,
        dtr_regions,
    );

    (start_areas, end_areas)
}

/// Filter positions and merge nearby ones into peak areas, with diagnostics.
///
/// This implements the same flow as filter_and_merge_to_areas but also
/// collects diagnostic information for each position that passes the first filter.
///
/// Returns (start_areas, end_areas, diagnostics) for classification and CSV export.
pub fn filter_and_merge_to_areas_with_diagnostics(
    reads_starts: &[u64],
    reads_ends: &[u64],
    left_clippings: &[u64],
    right_clippings: &[u64],
    coverage_reduced: &[u64],
    primary_reads: &[u64],
    unclipped_ratio: f64,
    _coverage_reduced_mean: f64,
    config: &PhageTerminiConfig,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (Vec<PeakArea>, Vec<PeakArea>, Vec<TerminusDiagnostic>) {
    let min_frequency = (config.min_frequency as f64) / 100.0;
    let min_events = config.min_events as u64;
    let max_distance_peaks = config.max_distance_peaks;
    let clipping_significance = config.clipping_significance;

    let mut diagnostics: Vec<TerminusDiagnostic> = Vec::new();

    // Step 1-2: Filter positions for starts (collecting diagnostics)
    let mut filtered_start_positions: Vec<(i32, u64, u64)> = Vec::new();
    for (idx, &spc) in reads_starts.iter().enumerate() {
        if spc == 0 {
            continue;
        }
        let pos = (idx + 1) as i32; // 1-indexed
        let cov = coverage_reduced.get(idx).copied().unwrap_or(0);
        // Sum clippings within max_distance_peaks of this position
        let clips = sum_clippings_in_range(left_clippings, idx, max_distance_peaks);
        let prim = primary_reads.get(idx).copied().unwrap_or(0);

        // Step 1: Local criteria (SPC >= 0.1*coverage_reduced AND SPC >= min_events)
        let local_threshold = (cov as f64 * min_frequency) as u64;
        if spc < local_threshold || spc < min_events {
            continue;  // Does NOT pass first filter, skip entirely
        }

        // Pre-filter clippings: only consider them significant if they meet threshold
        // This avoids false positives in low-coverage scenarios
        let effective_clips = if clips >= (prim as f64 * min_frequency) as u64 && clips >= min_events {
            clips
        } else {
            0  // Treat as no clippings - statistical test will auto-pass
        };

        // Position passes first filter - create diagnostic entry
        let observed_ratio = if spc + effective_clips > 0 {
            spc as f64 / (spc + effective_clips) as f64
        } else {
            1.0
        };
        let tau = if cov > 0 {
            spc as f64 / cov as f64
        } else {
            0.0
        };

        // Step 2: Statistical test (with reason tracking)
        let (passes_stat_test, stat_reason) = is_significant_terminus_with_reason(
            spc as u32, effective_clips as u32, unclipped_ratio, clipping_significance
        );

        let diag = TerminusDiagnostic {
            position: pos,
            terminus_type: "start".to_string(),
            spc,
            clippings: clips,  // Report raw clippings in diagnostic
            observed_ratio,
            primary_reads: prim,
            coverage_reduced: cov,
            tau,
            kept: passes_stat_test,  // Will be updated after global filter
            discarded_because: if passes_stat_test { String::new() } else { stat_reason },
        };

        if passes_stat_test {
            filtered_start_positions.push((pos, spc, effective_clips));
        }
        diagnostics.push(diag);
    }

    // Step 1-2: Filter positions for ends (collecting diagnostics)
    let mut filtered_end_positions: Vec<(i32, u64, u64)> = Vec::new();
    for (idx, &spc) in reads_ends.iter().enumerate() {
        if spc == 0 {
            continue;
        }
        let pos = (idx + 1) as i32;
        let cov = coverage_reduced.get(idx).copied().unwrap_or(0);
        // Sum clippings within max_distance_peaks of this position
        let clips = sum_clippings_in_range(right_clippings, idx, max_distance_peaks);
        let prim = primary_reads.get(idx).copied().unwrap_or(0);

        let local_threshold = (cov as f64 * 0.1) as u64;
        if spc < local_threshold || spc < min_events {
            continue;  // Does NOT pass first filter, skip entirely
        }

        // Pre-filter clippings: only consider them significant if they meet threshold
        let effective_clips = if clips >= (prim as f64 * 0.1) as u64 && clips >= min_events {
            clips
        } else {
            0  // Treat as no clippings - statistical test will auto-pass
        };

        // Position passes first filter - create diagnostic entry
        let observed_ratio = if spc + effective_clips > 0 {
            spc as f64 / (spc + effective_clips) as f64
        } else {
            1.0
        };
        let tau = if cov > 0 {
            spc as f64 / cov as f64
        } else {
            0.0
        };

        let (passes_stat_test, stat_reason) = is_significant_terminus_with_reason(
            spc as u32, effective_clips as u32, unclipped_ratio, clipping_significance
        );

        let diag = TerminusDiagnostic {
            position: pos,
            terminus_type: "end".to_string(),
            spc,
            clippings: clips,  // Report raw clippings in diagnostic
            observed_ratio,
            primary_reads: prim,
            coverage_reduced: cov,
            tau,
            kept: passes_stat_test,
            discarded_because: if passes_stat_test { String::new() } else { stat_reason },
        };

        if passes_stat_test {
            filtered_end_positions.push((pos, spc, effective_clips));
        }
        diagnostics.push(diag);
    }

    // Step 3: Merge nearby positions into areas (starts)
    let start_areas = merge_positions_to_areas(
        &filtered_start_positions,
        reads_starts,
        reads_ends,
        left_clippings,
        coverage_reduced,
        max_distance_peaks,
        genome_length,
        circular,
        dtr_regions,
    );

    // Step 3: Merge nearby positions into areas (ends)
    let end_areas = merge_positions_to_areas(
        &filtered_end_positions,
        reads_starts,
        reads_ends,
        right_clippings,
        coverage_reduced,
        max_distance_peaks,
        genome_length,
        circular,
        dtr_regions,
    );

    (start_areas, end_areas, diagnostics)
}

/// Merge filtered positions into peak areas based on proximity (DTR-aware).
fn merge_positions_to_areas(
    filtered_positions: &[(i32, u64, u64)], // (pos, spc, clips)
    reads_starts: &[u64],
    reads_ends: &[u64],
    _clippings: &[u64], // Not used directly - clips are in filtered_positions tuple
    coverage_reduced: &[u64],
    max_distance_peaks: i32,
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> Vec<PeakArea> {
    if filtered_positions.is_empty() {
        return Vec::new();
    }

    // Sort by position
    let mut positions: Vec<(i32, u64, u64)> = filtered_positions.to_vec();
    positions.sort_by_key(|(pos, _, _)| *pos);

    let mut areas: Vec<PeakArea> = Vec::new();

    // Start first area
    let (first_pos, first_spc, first_clips) = positions[0];
    let mut current_start = first_pos;
    let mut current_end = first_pos;
    let mut current_center = first_pos;
    let mut current_max_spc = first_spc;
    let mut current_total_spc = first_spc;
    let mut current_total_clips = first_clips;

    for &(pos, spc, clips) in positions.iter().skip(1) {
        let dist = dtr_aware_distance(
            current_end,
            pos,
            genome_length,
            circular,
            dtr_regions,
            max_distance_peaks,
        );

        if dist <= max_distance_peaks {
            // Extend current area
            current_end = pos;
            current_total_spc += spc;
            current_total_clips += clips;
            if spc > current_max_spc {
                current_max_spc = spc;
                current_center = pos;
            }
        } else {
            // Finalize current area and start new one
            let idx = (current_center - 1) as usize;
            let cov = coverage_reduced.get(idx).copied().unwrap_or(0) as u32;
            let starts = reads_starts.get(idx).copied().unwrap_or(0);
            let ends = reads_ends.get(idx).copied().unwrap_or(0);
            let tau = if cov > 0 {
                (starts + ends) as f64 / cov as f64
            } else {
                0.0
            };

            areas.push(PeakArea {
                start_pos: current_start,
                end_pos: current_end,
                center_pos: current_center,
                total_spc: current_total_spc as u32,
                total_clips: current_total_clips as u32,
                coverage: cov,
                tau,
            });

            // Start new area
            current_start = pos;
            current_end = pos;
            current_center = pos;
            current_max_spc = spc;
            current_total_spc = spc;
            current_total_clips = clips;
        }
    }

    // Finalize last area
    let idx = (current_center - 1) as usize;
    let cov = coverage_reduced.get(idx).copied().unwrap_or(0) as u32;
    let starts = reads_starts.get(idx).copied().unwrap_or(0);
    let ends = reads_ends.get(idx).copied().unwrap_or(0);
    let tau = if cov > 0 {
        (starts + ends) as f64 / cov as f64
    } else {
        0.0
    };

    areas.push(PeakArea {
        start_pos: current_start,
        end_pos: current_end,
        center_pos: current_center,
        total_spc: current_total_spc as u32,
        total_clips: current_total_clips as u32,
        coverage: cov,
        tau,
    });

    // Handle circular wrap-around merge
    if circular && areas.len() > 1 {
        // Copy values we need before mutating
        let first_end_pos = areas[0].end_pos;
        let first_total_spc = areas[0].total_spc;
        let first_total_clips = areas[0].total_clips;
        let first_center_pos = areas[0].center_pos;
        let first_coverage = areas[0].coverage;
        let first_tau = areas[0].tau;
        let first_start_pos = areas[0].start_pos;

        let last_idx = areas.len() - 1;
        let last_end_pos = areas[last_idx].end_pos;
        let last_total_spc = areas[last_idx].total_spc;
        let last_total_clips = areas[last_idx].total_clips;
        let last_center_pos = areas[last_idx].center_pos;
        let last_coverage = areas[last_idx].coverage;
        let last_tau = areas[last_idx].tau;
        let last_start_pos = areas[last_idx].start_pos;

        let wrap_dist = dtr_aware_distance(
            last_end_pos,
            first_start_pos,
            genome_length,
            circular,
            dtr_regions,
            max_distance_peaks,
        );

        if wrap_dist <= max_distance_peaks {
            // Merge first and last areas
            let merged_total_spc = first_total_spc + last_total_spc;
            let merged_total_clips = first_total_clips + last_total_clips;
            let (merged_center, merged_cov, merged_tau) = if first_total_spc >= last_total_spc {
                (first_center_pos, first_coverage, first_tau)
            } else {
                (last_center_pos, last_coverage, last_tau)
            };

            // Remove last area
            areas.pop();

            // Update first area to include merged data
            areas[0] = PeakArea {
                start_pos: last_start_pos, // Wrap around: start from last area
                end_pos: first_end_pos,
                center_pos: merged_center,
                total_spc: merged_total_spc,
                total_clips: merged_total_clips,
                coverage: merged_cov,
                tau: merged_tau,
            };
        }
    }

    areas
}

/// Check if a repeat is valid for phage termini analysis.
/// A valid repeat must have:
/// - One region starting within max_distance from position 1 (beginning)
/// - The other region ending within max_distance from contig_length (end)
pub fn is_valid_terminal_repeat(dup: &RepeatsData, contig_length: usize, max_distance: i32) -> bool {
    let first_start = dup.position1.min(dup.position2);
    let first_end = dup.position1.max(dup.position2);
    let second_start = dup.position1prime.min(dup.position2prime);
    let second_end = dup.position1prime.max(dup.position2prime);

    let contig_end = contig_length as i32;

    // Check if first region is at beginning and second region is at end
    let first_at_start = first_start <= max_distance;
    let second_at_end = second_end >= (contig_end - max_distance);

    // Check the reverse: first at end, second at beginning
    let first_at_end = first_end >= (contig_end - max_distance);
    let second_at_start = second_start <= max_distance;

    (first_at_start && second_at_end) || (first_at_end && second_at_start)
}


/// Translate a position from second DTR region to first DTR region.
/// Only applies to DTR (direct), not ITR (inverted).
/// Returns the translated position (or original if not in a DTR second region).
fn translate_to_first_dtr_region(pos: i32, dtr_regions: &[DtrRegion]) -> i32 {
    for dtr in dtr_regions {
        // Only translate for DTR (direct), not ITR (inverted)
        if !dtr.is_direct {
            continue;
        }
        // Check if pos is in second region
        if pos >= dtr.second_start && pos <= dtr.second_end {
            return dtr.first_start + pos - dtr.second_start;
        }
    }
    // Not in any DTR second region, return as-is
    pos
}

/// Deduplicate peaks that are DTR-equivalent and determine duplication status.
/// Groups peaks by equivalence and keeps the canonical (first region) position.
/// Only applies to DTR (direct repeats), not ITR (inverted repeats).
///
/// Returns (unique_peaks, duplication_status):
/// - unique_peaks: peaks with positions translated to first region
/// - duplication_status: Some(true) if all in DTR, Some(false) if all in ITR, None otherwise
fn deduplicate_dtr_equivalent(peaks: &[Peak], dtr_regions: &[DtrRegion]) -> (Vec<Peak>, Option<bool>) {
    if peaks.is_empty() {
        return (Vec::new(), None);
    }

    if dtr_regions.is_empty() {
        return (peaks.to_vec(), None);
    }

    let mut unique: Vec<Peak> = Vec::new();
    let mut all_in_dtr = true;
    let mut all_in_itr = true;
    let mut any_in_repeat = false;

    for peak in peaks {
        let pos = peak.position;

        // Check which region type this peak is in
        let mut in_dtr = false;
        let mut in_itr = false;

        for dtr in dtr_regions {
            let in_first = pos >= dtr.first_start && pos <= dtr.first_end;
            let in_second = pos >= dtr.second_start && pos <= dtr.second_end;

            if in_first || in_second {
                any_in_repeat = true;
                if dtr.is_direct {
                    in_dtr = true;
                } else {
                    in_itr = true;
                }
            }
        }

        if !in_dtr {
            all_in_dtr = false;
        }
        if !in_itr {
            all_in_itr = false;
        }

        // Translate to canonical (first region) position - only for DTR
        let canonical_pos = translate_to_first_dtr_region(pos, dtr_regions);

        // Check if this canonical position already exists
        let is_duplicate = unique.iter().any(|existing| {
            existing.position == canonical_pos
        });

        if !is_duplicate {
            unique.push(Peak {
                position: canonical_pos,
                value: peak.value,
            });
        }
    }

    // Determine duplication status
    let status = if !any_in_repeat {
        None
    } else if all_in_dtr && !all_in_itr {
        Some(true)  // All in DTR
    } else if all_in_itr && !all_in_dtr {
        Some(false) // All in ITR
    } else {
        None // Mixed or unclear
    };

    (unique, status)
}

/// Check if positions are in ITR regions: one in first, other in second.
/// Returns Some(repeat_size) if this is an ITR configuration, None otherwise.
fn check_itr_configuration(pos1: i32, pos2: i32, dtr_regions: &[DtrRegion]) -> Option<i32> {
    for dtr in dtr_regions {
        // Only check ITR regions (is_direct = false)
        if dtr.is_direct {
            continue;
        }

        let in_first_1 = pos1 >= dtr.first_start && pos1 <= dtr.first_end;
        let in_second_1 = pos1 >= dtr.second_start && pos1 <= dtr.second_end;
        let in_first_2 = pos2 >= dtr.first_start && pos2 <= dtr.first_end;
        let in_second_2 = pos2 >= dtr.second_start && pos2 <= dtr.second_end;

        // One position in first ITR region, other in second ITR region
        if (in_first_1 && in_second_2) || (in_second_1 && in_first_2) {
            let repeat_size = dtr.first_end - dtr.first_start + 1;
            return Some(repeat_size);
        }
    }
    None
}

/// Combine duplication status from start and end peaks.
/// Returns Some(true) if all in DTR, Some(false) if all in ITR, None otherwise.
fn combine_duplication_status(
    starts_status: Option<bool>,
    ends_status: Option<bool>,
    has_starts: bool,
    has_ends: bool,
) -> Option<bool> {
    match (starts_status, ends_status, has_starts, has_ends) {
        // Both have same status
        (Some(true), Some(true), _, _) => Some(true),   // All DTR
        (Some(false), Some(false), _, _) => Some(false), // All ITR
        // Only one side has peaks - use that status
        (s, _, true, false) => s,
        (_, s, false, true) => s,
        // Mixed or no peaks
        _ => None,
    }
}

/// Classify phage packaging mechanism based on peak configuration.
/// Returns (mechanism, left_termini, right_termini, duplication_status).
///
/// DTR equivalence is handled here: peaks from both DTR regions are received,
/// and we deduplicate equivalent peaks to determine the "real" peak count
/// for classification. Termini positions are returned as-is (no expansion).
pub fn classify_packaging(
    start_peaks: &[Peak],
    end_peaks: &[Peak],
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (String, Vec<i32>, Vec<i32>, Option<bool>) {
    // Deduplicate DTR-equivalent peaks and get duplication status
    let (unique_starts, starts_status) = deduplicate_dtr_equivalent(start_peaks, dtr_regions);
    let (unique_ends, ends_status) = deduplicate_dtr_equivalent(end_peaks, dtr_regions);

    // Use unique counts for classification logic
    match (unique_starts.len(), unique_ends.len()) {
        (0, 0) => ("No_packaging".to_string(), vec![], vec![], None),

        (1, 0) => {
            // Only start peak - PAC
            // Return all original positions (not deduplicated)
            let left: Vec<i32> = start_peaks.iter().map(|p| p.position).collect();
            ("PAC".to_string(), left, vec![], starts_status)
        }

        (0, 1) => {
            // Only end peak - PAC
            let right: Vec<i32> = end_peaks.iter().map(|p| p.position).collect();
            ("PAC".to_string(), vec![], right, ends_status)
        }

        (1, 1) => {
            // Two unique peaks - classify based on distance and order
            // Use the first unique peak positions for classification
            let start_pos = unique_starts[0].position;
            let end_pos = unique_ends[0].position;

            // Calculate actual genomic distance (shortest path in circular)
            let distance = if circular {
                circular_distance(start_pos, end_pos, genome_length)
            } else {
                (start_pos - end_pos).abs()
            };

            let genome_10pct = (genome_length as f64 * 0.1) as i32;

            // Determine order: which comes first going forward (clockwise in circular)
            // For linear: simply check positions
            // For circular: check which direction is shorter to go from start to end
            let end_before_start = if circular {
                // Calculate forward distance from start to end
                let gl = genome_length as i32;
                let forward_dist = if end_pos >= start_pos {
                    end_pos - start_pos
                } else {
                    (gl - start_pos) + end_pos
                };
                // If forward distance equals circular_distance, end is ahead
                // If they're different, end is behind (going backward is shorter)
                forward_dist != distance
            } else {
                end_pos < start_pos
            };

            // Check for ITR configuration: COS-like peaks in ITR regions
            // For ITR, we use the repeat size from autoblast, not the peak distance
            let itr_config = if distance <= 20 {
                check_itr_configuration(start_pos, end_pos, dtr_regions)
            } else {
                None
            };

            // Determine suffix and position ordering based on directionality
            let suffix = if end_before_start { "_3'" } else { "_5'" };
            let (_left_pos, _right_pos) = if end_before_start {
                (end_pos, start_pos)
            } else {
                (start_pos, end_pos)
            };

            let mechanism = if let Some(repeat_size) = itr_config {
                // ITR detected - classify based on repeat size from autoblast
                if repeat_size <= 1000 {
                    format!("ITR_short{}", suffix)
                } else if repeat_size <= genome_10pct {
                    format!("ITR_long{}", suffix)
                } else {
                    format!("ITR_outlier{}", suffix)
                }
            } else if distance < 2 {
                // Very close - cohesive ends (no suffix)
                "COS".to_string()
            } else if distance <= 20 {
                // Close cohesive ends with directionality
                format!("COS{}", suffix)
            } else if distance <= 1000 {
                // Short direct terminal repeats
                format!("DTR_short{}", suffix)
            } else if distance <= genome_10pct {
                // Long direct terminal repeats
                format!("DTR_long{}", suffix)
            } else {
                // Very distant - outlier
                format!("DTR_outlier{}", suffix)
            };

            // Return all original positions (not just deduplicated ones)
            let mut left: Vec<i32> = start_peaks.iter().map(|p| p.position).collect();
            let mut right: Vec<i32> = end_peaks.iter().map(|p| p.position).collect();
            left.sort();
            right.sort();
            let dup_status = combine_duplication_status(starts_status, ends_status, true, true);
            (mechanism, left, right, dup_status)
        }

        _ => {
            // Multiple unique peaks - Unknown_packaging
            // Return all original positions
            let mut left_termini: Vec<i32> = start_peaks.iter().map(|p| p.position).collect();
            let mut right_termini: Vec<i32> = end_peaks.iter().map(|p| p.position).collect();

            // Sort for consistent output
            left_termini.sort();
            right_termini.sort();

            let dup_status = combine_duplication_status(starts_status, ends_status, !unique_starts.is_empty(), !unique_ends.is_empty());
            ("Unknown_packaging".to_string(), left_termini, right_termini, dup_status)
        }
    }
}

/// Deduplicate peak areas that are DTR-equivalent and determine duplication status.
/// Groups areas by equivalence (using center_pos) and keeps the canonical (first region) position.
/// Only applies to DTR (direct repeats), not ITR (inverted repeats).
///
/// Returns (unique_areas, duplication_status):
/// - unique_areas: areas with center_pos translated to first region
/// - duplication_status: Some(true) if all in DTR, Some(false) if all in ITR, None otherwise
fn deduplicate_dtr_equivalent_areas(areas: &[PeakArea], dtr_regions: &[DtrRegion]) -> (Vec<PeakArea>, Option<bool>) {
    if areas.is_empty() {
        return (Vec::new(), None);
    }

    if dtr_regions.is_empty() {
        return (areas.to_vec(), None);
    }

    let mut unique: Vec<PeakArea> = Vec::new();
    let mut all_in_dtr = true;
    let mut all_in_itr = true;
    let mut any_in_repeat = false;

    for area in areas {
        let pos = area.center_pos;

        // Check which region type this area's center is in
        let mut in_dtr = false;
        let mut in_itr = false;

        for dtr in dtr_regions {
            let in_first = pos >= dtr.first_start && pos <= dtr.first_end;
            let in_second = pos >= dtr.second_start && pos <= dtr.second_end;

            if in_first || in_second {
                any_in_repeat = true;
                if dtr.is_direct {
                    in_dtr = true;
                } else {
                    in_itr = true;
                }
            }
        }

        if !in_dtr {
            all_in_dtr = false;
        }
        if !in_itr {
            all_in_itr = false;
        }

        // Translate center_pos to canonical (first region) position - only for DTR
        let canonical_pos = translate_to_first_dtr_region(pos, dtr_regions);

        // Check if this canonical position already exists (within tolerance)
        let is_duplicate = unique.iter().any(|existing| {
            (existing.center_pos - canonical_pos).abs() <= 20 // Use 20bp tolerance
        });

        if !is_duplicate {
            // Store area with canonical position for correct distance calculation
            let mut translated_area = area.clone();
            translated_area.center_pos = canonical_pos;
            unique.push(translated_area);
        }
    }

    // Determine duplication status
    let status = if !any_in_repeat {
        None
    } else if all_in_dtr && !all_in_itr {
        Some(true)  // All in DTR
    } else if all_in_itr && !all_in_dtr {
        Some(false) // All in ITR
    } else {
        None // Mixed or unclear
    };

    (unique, status)
}

/// Convert PeakArea to TerminusArea for database storage.
fn peak_area_to_terminus_area(area: &PeakArea) -> TerminusArea {
    TerminusArea {
        start_pos: area.start_pos,
        end_pos: area.end_pos,
        center_pos: area.center_pos,
        total_spc: area.total_spc,
        coverage: area.coverage,
        tau: area.tau,
    }
}

/// Classify phage packaging mechanism based on peak areas.
/// Returns (mechanism, left_termini, right_termini, duplication_status).
///
/// This is the new area-based version of classify_packaging.
/// DTR equivalence is handled here: areas from both DTR regions are received,
/// and we deduplicate equivalent areas to determine the "real" area count
/// for classification.
pub fn classify_packaging_areas(
    start_areas: &[PeakArea],
    end_areas: &[PeakArea],
    genome_length: usize,
    circular: bool,
    dtr_regions: &[DtrRegion],
) -> (String, Vec<TerminusArea>, Vec<TerminusArea>, Option<bool>) {
    // Deduplicate DTR-equivalent areas and get duplication status
    let (unique_starts, starts_status) = deduplicate_dtr_equivalent_areas(start_areas, dtr_regions);
    let (unique_ends, ends_status) = deduplicate_dtr_equivalent_areas(end_areas, dtr_regions);

    // Use unique counts for classification logic
    match (unique_starts.len(), unique_ends.len()) {
        (0, 0) => ("No_packaging".to_string(), vec![], vec![], None),

        (1, 0) => {
            // Only start area - PAC
            let left: Vec<TerminusArea> = start_areas.iter().map(peak_area_to_terminus_area).collect();
            ("PAC".to_string(), left, vec![], starts_status)
        }

        (0, 1) => {
            // Only end area - PAC
            let right: Vec<TerminusArea> = end_areas.iter().map(peak_area_to_terminus_area).collect();
            ("PAC".to_string(), vec![], right, ends_status)
        }

        (1, 1) => {
            // Two unique areas - classify based on distance and order
            // Use center positions for classification
            let start_pos = unique_starts[0].center_pos;
            let end_pos = unique_ends[0].center_pos;

            // Calculate actual genomic distance (shortest path in circular)
            let distance = if circular {
                circular_distance(start_pos, end_pos, genome_length)
            } else {
                (start_pos - end_pos).abs()
            };

            let genome_10pct = (genome_length as f64 * 0.1) as i32;

            // Determine order: which comes first going forward (clockwise in circular)
            let end_before_start = if circular {
                let gl = genome_length as i32;
                let forward_dist = if end_pos >= start_pos {
                    end_pos - start_pos
                } else {
                    (gl - start_pos) + end_pos
                };
                forward_dist != distance
            } else {
                end_pos < start_pos
            };

            // Check for ITR configuration
            let itr_config = if distance <= 20 {
                check_itr_configuration(start_pos, end_pos, dtr_regions)
            } else {
                None
            };

            let suffix = if end_before_start { "_3'" } else { "_5'" };

            let mechanism = if let Some(repeat_size) = itr_config {
                if repeat_size <= 1000 {
                    format!("ITR_short{}", suffix)
                } else if repeat_size <= genome_10pct {
                    format!("ITR_long{}", suffix)
                } else {
                    format!("ITR_outlier{}", suffix)
                }
            } else if distance < 2 {
                "COS".to_string()
            } else if distance <= 20 {
                format!("COS{}", suffix)
            } else if distance <= 1000 {
                format!("DTR_short{}", suffix)
            } else if distance <= genome_10pct {
                format!("DTR_long{}", suffix)
            } else {
                format!("DTR_outlier{}", suffix)
            };

            // Return all original areas (not just deduplicated ones)
            let left: Vec<TerminusArea> = start_areas.iter().map(peak_area_to_terminus_area).collect();
            let right: Vec<TerminusArea> = end_areas.iter().map(peak_area_to_terminus_area).collect();
            let dup_status = combine_duplication_status(starts_status, ends_status, true, true);
            (mechanism, left, right, dup_status)
        }

        _ => {
            // Multiple unique areas - Unknown_packaging
            let left_termini: Vec<TerminusArea> = start_areas.iter().map(peak_area_to_terminus_area).collect();
            let right_termini: Vec<TerminusArea> = end_areas.iter().map(peak_area_to_terminus_area).collect();

            let dup_status = combine_duplication_status(starts_status, ends_status, !unique_starts.is_empty(), !unique_ends.is_empty());
            ("Unknown_packaging".to_string(), left_termini, right_termini, dup_status)
        }
    }
}

