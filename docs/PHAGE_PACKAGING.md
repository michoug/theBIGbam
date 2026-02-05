# Phage Packaging Classification Algorithm

The `processing_phage_packaging.rs` module classifies phage packaging mechanisms through a multi-step filtering and classification pipeline.

## 0. Computing Read Starts and Ends (SPC)

The phage termini detection relies on **reads_starts** and **reads_ends** arrays that count where reads begin and end along the contig.

### Filtering: Only Matching Termini

A read's terminus is only counted if it starts with a **match** (no soft-clipping or mismatch). This is verified using both the CIGAR string and MD tag. Reads starting with clipping or mismatches are excluded because they likely represent sequencing artifacts or misalignments rather than true biological termini.

### Short vs Long Reads

**Short reads** (paired-end): The 5' end of read1 (+ strand) represents the fragment start; the 5' end of read2 (- strand) represents the fragment end. Thus:

- reads_starts = start positions of + strand reads
- reads_ends = start positions of - strand reads (which are fragment 3' ends)

**Long reads**: Each read is **split in half** and each terminus is evaluated independently. This prevents losing ~80% of reads that have clipping at one end but a clean match at the other. Both strands contribute to both arrays:

- reads_starts = 5' ends from both strands (start_plus + start_minus)
- reads_ends = 3' ends from both strands (end_plus + end_minus)

### SPC: Starts Per Coverage

**SPC (Starts Per Coverage)** = reads_starts / coverage_reduced at each position. This normalizes the read start count by coverage, allowing detection of positions where reads start disproportionately often. Similarly, **EPC (Ends Per Coverage)** = reads_ends / coverage_reduced.

A high SPC/EPC indicates a potential packaging terminus: reads preferentially start/end at that position rather than being uniformly distributed.

## 1. Contig Selection

Only contigs where ≥ **min_aligned_fraction** (default 90%) of positions have at least one aligned read are processed.

## 3. Position Filtering

Each position with read_starts > 0 or read_ends > 0 is evaluated:

Position must pass BOTH criteria:

- SPC ≥ **min_frequency** × coverage_reduced (default 10%)
- SPC ≥ **min_events** (default 10)

Where SPC = read_starts (for start positions) or read_ends (for end positions).

## 4. Peak Merging

Positions passing the filter are merged into peak areas:

- Nearby positions within **max_distance_peaks** (default 20bp) merge together
- Distance calculations are **DTR/ITR-aware**: positions in second repeat region are treated as equivalent to their first region counterparts
- For circular genomes, wrap-around merging occurs
- **Positions are stored as canonical** (translated to first DTR region for DTR; ITR positions remain as-is)
- Each area stores: start_pos, end_pos, center_pos (highest SPC), total_spc
2. Terminal Repeat Detection (from autoblast)

Terminal repeats (DTR/ITR) are identified from autoblast results using:

- Identity ≥ **min_identity_dtr** (default 90%)
- One region within **max_distance_duplication** (default 100bp) of contig start
- Other region within **max_distance_duplication** of contig end

These regions enable:

- **DTR-aware distance calculations** for peak merging and clipping aggregation
- **Canonical position storage** (peaks translated to first region)
- **Both-copies confirmation** for clipping validation (DTR only)

## 5. Clipping Pre-filter

Individual clipping positions are pre-filtered before aggregation. A clipping position is significant if:

- clippings ≥ **min_frequency** × primary_reads (default 10%)
- clippings ≥ **min_events** (default 10)

For start areas: evaluate left_clippings at each position  
For end areas: evaluate right_clippings at each position

### 5b. DTR Both-Copies Confirmation

For **DTR only** (not ITR): real biological clippings appear at both DTR copies at equivalent positions. Boundary artifacts (reads can't extend past contig edge) appear at only one copy. 

Any clipping that lacks a counterpart at the equivalent position in the other DTR copy is zeroed out. This prevents false positives from contig boundary effects.

## 6. Area Clipping Aggregation

For each peak area, sum all pre-filtered clippings:

- Within the area bounds (start_pos to end_pos)
- OR within **max_distance_peaks** of area edges

**Distance awareness differs by repeat type:**

| Repeat Type | Distance Calculation                 | Rationale                                                                            |
| ----------- | ------------------------------------ | ------------------------------------------------------------------------------------ |
| **DTR**     | DTR-aware + circular-aware           | Clippings at equivalent positions in both copies represent the same signal           |
| **ITR**     | Circular-aware only (no ITR mapping) | Left clippings in first ITR ≠ left clippings in second ITR (orientation is inverted) |
| **None**    | Circular-aware only                  | Standard distance calculation                                                        |

## 7. Statistical Test

Tests whether a peak area has **significantly more clippings** than expected globally.

**Expected clippings calculation:**

- `clipped_ratio = (primary_reads_mean - coverage_reduced_mean) / coverage_reduced_mean`
- `expected_clippings = clipped_ratio × total_SPC`

**Decision logic:**

- If `sum_clippings > expected_clippings` significantly → **DISCARD** (excess clippings indicate artifact)
- Otherwise → **KEEP**

Statistical significance is determined using **clipping_significance** (default 0.05 → z_critical = 1.645).

## 8. DTR Deduplication

Before classification, remaining peak areas are deduplicated:

- Areas already store canonical positions (from step 4)
- Duplicate areas at the same canonical position (within 20bp tolerance) are removed
- Result: unique peak count for classification

**Important:** Canonical positions ensure correct distance calculation for classification (e.g., T7 phage: start at pos 1, end at pos 160 → DTR distance = 159bp, not ~39000bp if using second copy position).

## 9. Classification

Based on unique peak count after deduplication:

| Start Areas | End Areas | Classification                 |
| ----------- | --------- | ------------------------------ |
| 0           | 0         | No_packaging                   |
| 1           | 0         | PAC                            |
| 0           | 1         | PAC                            |
| 1           | 1         | See distance-based rules below |
| >1 or >1    |           | Unknown_packaging              |

### Distance-based Classification (1 start, 1 end)

For ITR configuration (both peaks in ITR regions, distance ≤20bp):

- ITR length ≤1000bp → **ITR_short_5'/3'**
- ITR length ≤10% genome → **ITR_long_5'/3'**
- ITR length >10% genome → **ITR_outlier_5'/3'**

Otherwise, based on distance between canonical positions:

- <2bp → **COS** (blunt cohesive ends)
- ≤20bp → **COS_5'/3'** (cohesive with overhang)
- ≤1000bp → **DTR_short_5'/3'**
- ≤10% genome → **DTR_long_5'/3'**
- > 10% genome → **DTR_outlier_5'/3'**

The suffix (_5' or _3') indicates overhang orientation based on whether end comes before or after start in genomic coordinates.

## 10. Output

- **Mechanism**: Classification string (e.g., "DTR_short_5'")
- **Left termini**: All start area positions (original positions, including both DTR copies if present)
- **Right termini**: All end area positions (original positions, including both DTR copies if present)
- **Duplication status**: `true` if all peaks in DTR, `false` if all in ITR, `null` otherwise
- **Repeat length**: Distance between start and end peak centers (when exactly 1 unique of each after deduplication)

## 11. Diagnostic CSV Output

A diagnostic CSV file (`{db_name}_phagetermini.csv`) is generated with one row per peak area. Columns:

| Column                | Description                                                                                |
| --------------------- | ------------------------------------------------------------------------------------------ |
| Contig_name           | Contig name                                                                                |
| Sample_name           | Sample name                                                                                |
| Clipped_ratio         | Global clipped ratio: (primary_reads_mean - coverage_reduced_mean) / coverage_reduced_mean |
| Coverage_reduced_mean | Global mean coverage_reduced                                                               |
| Type                  | "start" or "end"                                                                           |
| Position              | Center position of peak area (1-indexed, original position)                                |
| SPC                   | Total SPC in peak area                                                                     |
| Clippings             | Sum of pre-filtered clippings in/near area                                                 |
| Expected_clippings    | clipped_ratio × total_SPC                                                                  |
| Kept                  | "yes" or "no"                                                                              |
| Discarded_because     | Reason if discarded, empty if kept                                                         |

Discard reasons:

- `excess_clippings: z=X.XX > z_critical=X.XX` - significantly more clippings than expected

## Summary: DTR vs ITR Handling

| Aspect                   | DTR (Direct)                  | ITR (Inverted)                               |
| ------------------------ | ----------------------------- | -------------------------------------------- |
| Peak merging             | Simple distance               | Simple distance                              |
| Position storage         | Original positions            | Original positions                           |
| Both-copies confirmation | Yes (zeros artifacts)         | No                                           |
| Clipping aggregation     | Simple distance               | Simple distance                              |
| Deduplication            | At classification (canonical) | At classification (no merge)                 |
| ITR repeat length        | N/A                           | From autoblast (first_end - first_start + 1) |

**Simplified architecture**: Peak areas are processed independently until classification. DTR-equivalent areas (from both copies of the repeat) are stored separately in the database, then deduplicated at classification time using canonical positions. This ensures all original peak information is preserved while correctly classifying DTR/ITR mechanisms.

The key difference at classification: DTR copies are directly equivalent (same sequence, same orientation), so areas from both regions are merged to a single canonical area for counting. ITR copies are inverted, so one peak in each region represents a legitimate COS-like configuration, not duplication.
