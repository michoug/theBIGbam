# Phage Packaging Classification Algorithm

The `processing_phage_packaging.rs` module classifies phage packaging mechanisms through a multi-step filtering and classification pipeline.

## 1. Computing Read Starts and Ends

The phage termini detection relies on **reads_starts** and **reads_ends** arrays that count where reads begin and end along the contig.

### Filtering: Only Matching Termini

A read's terminus is only counted if it starts with a **match** (no soft-clipping or mismatch). This is verified using both the CIGAR string and MD tag. Reads starting with clipping or mismatches are excluded because they likely represent sequencing artifacts or misalignments rather than true biological termini.

### Short vs Long Reads

**Short reads**: For paired-end sequencing, the 5' end of read1 (+ strand) represents the fragment start; the 5' end of read2 (- strand) represents the fragment end. Thus:

- reads_starts = start positions of + strand reads
- reads_ends = start positions of - strand reads (which are fragment 3' ends)

For single-end sequencing, only one fragment terminus is observed. The read is classified as `reads_starts` if it maps to the + strand and as `reads_ends` if it maps to the − strand.

**Long reads**: Each read is **split in half** and each terminus is evaluated independently. This prevents losing ~80% of reads that have clipping at one end but a clean match at the other. Both strands contribute to both arrays:

- reads_starts = 5' ends from both strands (start_plus + start_minus)
- reads_ends = 3' ends from both strands (end_plus + end_minus)

## 2. Contig Selection

Only contigs where ≥ **min_aligned_fraction** (default 90%) of positions have at least one aligned read are processed.

## 3. Terminal repeats (DTR/ITR) are identified:

Autoblast results are filtered to keep real terminal repeats:

- Identity ≥ **min_identity_dtr** (default 90%)
- One region within **max_distance_duplication** (default 100bp) of contig start
- Other region within **max_distance_duplication** of contig end

These regions enable:

- **Both-copies confirmation** for clipping validation (DTR only)
- Deduplication of peaks at equivalent positions in both DTR copies (at classification time)

## 3. Position Filtering

A position must pass two criteria to be considered as potential termini:

- read_starts ≥ **min_frequency** × coverage_reduced (default 10%)
- read_starts ≥ **min_events** (default 10)

Where the raw count used is read_starts (for start positions) or read_ends (for end positions).

## 4. Peak Merging

Positions passing the filter are merged into peak areas:

- Nearby positions within **max_distance_peaks** (default 20bp) merge together
- For circular genomes, wrap-around merging occurs
- Each area stores: start_pos, end_pos, peak_size, center_pos (position with highest signal), total_spc (spc stands for Starting Position Coverage)

## 5. Peak Testing (Poisson Significance)

After merging, each peak area is tested for statistical significance using a Poisson model. Under uniform random fragmentation, read starts at each position follow a Poisson distribution with rate proportional to local coverage.

**Background rate computation:**

- `pstart = sum(reads_starts) / sum(coverage_reduced)` (computed separately for starts and ends)

**Per-window expected count:**

- The window spans from `start_pos` to `end_pos` of the merged area
- `λ_W = sum(coverage_reduced[i] × pstart)` for all positions i in the window
- This accounts for variable coverage across the window

**Size:**

- Each area's `size` field is computed as the distance between `start_pos` and `end_pos` + 1

**Observed count:**

- `X_W = area.total_spc` (total read starts/ends in the window)

**Significance test:**

- Model: `X_W ~ Poisson(λ_W)`
- p-value: `P(X >= X_W | Poisson(λ_W))` (upper-tail, via `statrs` crate)
- Bonferroni correction: `adjusted_pvalue = pvalue × K` where K = number of windows tested (separate K for starts and ends)
- Areas with `adjusted_pvalue > 0.05` have `passed_poisson_test` set to false

**Wrap-around handling:** For circular genomes where `start_pos > end_pos`, the expected count sums coverage from `start_pos` to genome end plus position 1 to `end_pos`.

The Poisson test validates that the peak itself is significant (more starts than expected under uniform fragmentation). The subsequent clipping test checks whether a significant peak is an artifact (excess clippings indicating misalignment). These are complementary: the Poisson test filters noise, the clipping test filters artifacts. An area must pass **both** tests (`passed_poisson_test && passed_clipping_test`) to be kept.

## 6. Clipping Testing

### 6a. Clipping Pre-filter

Individual clipping positions are pre-filtered before aggregation. A clipping position is significant if:

- clippings ≥ **min_frequency** × primary_reads (default 10%)
- clippings ≥ **min_events** (default 10)

For start areas: evaluate left_clippings at each position
For end areas: evaluate right_clippings at each position

### 6b. DTR Both-Copies Confirmation

For **DTR only** (not ITR): real biological clippings appear at both DTR copies at equivalent positions. Boundary artifacts (reads can't extend past contig edge) appear at only one copy.

Any clipping that lacks a counterpart at the equivalent position in the other DTR copy is zeroed out. This prevents false positives from contig boundary effects.

### 6c. Area Clipping Aggregation

For each peak area, sum all pre-filtered clippings:

- Within the area bounds (start_pos to end_pos)
- OR within **max_distance_peaks** of area edges

**Distance calculation uses circular-aware distance only** (not DTR-aware). DTR deduplication happens at classification time, not during clipping aggregation.

| Repeat Type | Distance Calculation | Rationale                                                                            |
| ----------- | -------------------- | ------------------------------------------------------------------------------------ |
| **DTR**     | Circular-aware only  | DTR deduplication handled separately at classification                               |
| **ITR**     | Circular-aware only  | Left clippings in first ITR ≠ left clippings in second ITR (orientation is inverted) |
| **None**    | Circular-aware only  | Standard distance calculation                                                        |

### 6d. Statistical Test

Tests whether a peak area has **significantly more clippings** than expected globally.

**Expected clippings calculation:**

- `clipped_ratio = (primary_reads_mean - coverage_reduced_mean) / coverage_reduced_mean`
- `expected_clippings = clipped_ratio × total_spc`

**Decision logic:**

- If `sum_clippings > expected_clippings` significantly → **DISCARD** (`passed_clipping_test = false`)
- Otherwise → **KEEP**

Statistical significance is determined using **clipping_significance** (default 0.05 → z_critical = 1.645).

## 7. DTR Deduplication

Before classification, remaining peak areas are deduplicated:

- Each area's `center_pos` is translated to the first DTR region (canonical position)
- Duplicate areas at the same canonical position (within 20bp tolerance) are removed
- Result: unique peak count for classification

**Important:** Canonical positions ensure correct distance calculation for classification (e.g., T7 phage: start at pos 1, end at pos 160 → DTR distance = 159bp, not ~39000bp if using second copy position).

## 8. Classification

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

## 9. Output

- **Mechanism**: Classification string (e.g., "DTR_short_5'")
- **Left termini**: All start area positions (original positions)
- **Right termini**: All end area positions (original positions)
- **Duplication status**: `true` if all peaks in DTR, `false` if all in ITR, `null` otherwise
- **Repeat length**: Distance between start and end peak centers (when exactly 1 unique of each after deduplication)

## Summary: DTR vs ITR Handling

| Aspect                   | DTR (Direct)                  | ITR (Inverted)                               |
| ------------------------ | ----------------------------- | -------------------------------------------- |
| Peak merging             | Simple distance               | Simple distance                              |
| Position storage         | Original positions            | Original positions                           |
| Both-copies confirmation | Yes (zeros artifacts)         | No                                           |
| Clipping aggregation     | Circular-aware only           | Circular-aware only                          |
| Deduplication            | At classification (canonical) | At classification (no merge)                 |
| ITR repeat length        | N/A                           | From autoblast (first_end - first_start + 1) |

**Simplified architecture**: Peak areas are processed independently until classification. DTR-equivalent areas (from both copies of the repeat) are stored separately in the database, then deduplicated at classification time using canonical positions. This ensures all original peak information is preserved while correctly classifying DTR/ITR mechanisms.

The key difference at classification: DTR copies are directly equivalent (same sequence, same orientation), so areas from both regions are merged to a single canonical area for counting. ITR copies are inverted, so one peak in each region represents a legitimate COS-like configuration, not duplication.
