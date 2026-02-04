# Phage Packaging Classification Algorithm

The processing_phage_packaging.rs file classifies phage packaging mechanisms
through a multi-step pipeline:

## Contig selection:

Only contigs/sample pairs where more than **min_aligned_fraction** bp are covered
by at least one read are considered (default 90%). This is calculated on the
original coverage data before any processing.

*processing.rs: lines 537-542*

## Terminal repeat detection:

Terminal repeats (DTR/ITR) are identified using 2 criteria:

- Both regions share ≥ **min_identity_dtr** identity (default 90%)

- Both regions are located close to each contig end, at less than **max_distance_duplication** (default 100bp)

These regions are used for DTR-aware distance calculations and peak deduplication,
but the original signal data is preserved (no merging or zeroing).

*processing.rs: lines 247-282*

## Peak filtering (filter_peaks_by_clippings function):

Read start and end peaks are filtered based on three criteria:

- They must be more frequent than **min_frequency** compared to local coverage
(default 0.1)

- They must occur more than **min_events** times (default 10)

- They must NOT have nearby clipping events within **max_distance_peaks** (default 20bp). Clipping events indicate alignment artifacts rather than true termini,
so peaks near clippings are discarded

Distance calculations are **DTR/ITR-aware**: a clipping in the second repeat region
is considered "near" a peak at the equivalent position in the first region.
For DTR (direct): `second_start + offset` maps to `first_start + offset`.
For ITR (inverted): `second_start + offset` maps to `first_end - offset`.

*processing.rs: lines 553-563*

## Peak merging (merge_peaks function):

Nearby peaks within **max_distance_peaks** (default 20bp) are merged into a single
peak, keeping the peak with the strongest signal. Distance calculations use
**dtr_aware_distance**, which treats positions in the second repeat region as
equivalent to their counterparts in the first region (with proper DTR vs ITR mapping).
For circular genomes, wrap-around merging also occurs.

*processing.rs: lines 565-578*

## Packaging classification (classify_packaging function):

Classification proceeds in two steps:

### Step 1: Deduplicate DTR-equivalent peaks

The **deduplicate_dtr_equivalent** function:
- Translates peaks from second DTR region to canonical (first region) positions
- Removes duplicate peaks at the same canonical position
- Tracks whether peaks are in DTR regions, ITR regions, or neither
- Only applies to DTR (direct repeats); ITR peaks are not deduplicated

This determines the "real" unique peak count for classification.

### Step 2: Classify based on unique peak count

Based on the number and configuration of unique peaks:

### (0 peaks in reads_starts, 0 peaks in read_ends) = "No_packaging"

### (1, 0) or (0, 1) = **"PAC"** (headful packaging with single terminus)

### (1, 1) = Classification depends on configuration and distance

If both termini are part of an ITR (**check_itr_configuration function**), they
appear very close to each other on the sequence. In this case, the ITR length
from autoblast is used instead of the physical distance between peaks:

- ≤1000bp → **"ITR_short_5'/3'"** (short inverted terminal repeats)

- ≤10% genome → **"ITR_long_5'/3'"** (long inverted terminal repeats)

- >10% genome → **"ITR_outlier_5'/3'"**

Otherwise, the distance between peaks is used:

- <2bp → **"COS"** (cohesive ends, blunt)

- ≤20bp → **"COS_5'" or "COS_3'"** (cohesive with overhang directionality)

- ≤1000bp → **"DTR_short_5'/3'"** (short direct terminal repeats)

- ≤10% genome → **"DTR_long_5'/3'"** (long direct terminal repeats)

- >10% genome → **"DTR_outlier_5'/3'"**

### Other configurations → **"Unknown_packaging"**

## Output:

For all configurations:

- **Termini positions**: All original peak positions are returned (including both
  first and second DTR region positions if peaks exist in both). No expansion or
  translation is applied to the output.

- **Duplication status**: Determined during deduplication by tracking which regions
  peaks belong to. The **combine_duplication_status** function merges the status
  from start and end peaks: `Some(true)` if all peaks are in DTR regions,
  `Some(false)` if all are in ITR regions, `None` otherwise.

- **Suffix (_5' or _3')**: Indicates overhang orientation, determined by whether
  the end peak comes before or after the start peak in genomic coordinates
  (accounting for circularity).
