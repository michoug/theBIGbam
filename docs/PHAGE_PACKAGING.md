# Phage Packaging Classification Algorithm

The `processing_phage_packaging.rs` module classifies phage packaging mechanisms through a multi-step filtering and classification pipeline.

## 1. Contig Selection

Only contigs where ≥ **min_aligned_fraction** (default 90%) of positions have at least one aligned read are processed.

## 2. Terminal Repeat Detection (from autoblast)

Terminal repeats (DTR/ITR) are identified from autoblast results using:
- Identity ≥ **min_identity_dtr** (default 90%)
- One region within **max_distance_duplication** (default 100bp) of contig start
- Other region within **max_distance_duplication** of contig end

These regions enable DTR-aware distance calculations and peak deduplication.

## 3. Position Filtering

Each position with read_starts > 0 or read_ends > 0 is evaluated:

### Step 3a: Local Filter
Position must pass BOTH criteria:
- SPC ≥ **min_frequency** × coverage_reduced (default 10%)
- SPC ≥ **min_events** (default 10)

Where SPC = read_starts (for start positions) or read_ends (for end positions).

### Step 3b: Clipping Aggregation
Clippings are summed within **max_distance_peaks** (default 20bp) of the peak position. This ensures that nearby clipping events are properly associated with the peak they relate to.

For start positions: sum left_clippings within ±max_distance_peaks  
For end positions: sum right_clippings within ±max_distance_peaks

### Step 3c: Clipping Pre-filter
Aggregated clippings are only considered significant if they meet threshold:
- clippings ≥ **min_frequency** × primary_reads (default 10%)
- clippings ≥ **min_events** (default 10)

If clippings don't meet threshold, they're treated as 0 (position auto-passes statistical test).

### Step 3d: Statistical Test
Tests whether position has **significantly more clippings** than expected globally.

- `observed_ratio = SPC / (SPC + effective_clippings)`
- `expected_ratio = unclipped_ratio` (global proportion of unclipped reads)

**Decision logic:**
- If `observed_ratio ≥ expected_ratio` → **KEEP** (same or fewer clippings than average)
- If `observed_ratio < expected_ratio` → test if significantly worse:
  - z = (expected - observed) / std_error
  - If z > z_critical → **DISCARD** (excess clippings indicate artifact)
  - Otherwise → **KEEP** (not significantly worse)

The z_critical depends on **clipping_significance** (default 0.05 → z_critical = 1.645).

## 4. Peak Merging

Positions passing all filters are merged into peak areas:
- Nearby positions within **max_distance_peaks** (default 20bp) merge together
- Distance calculations are **DTR-aware**: positions in second repeat region are treated as equivalent to their first region counterparts
- For circular genomes, wrap-around merging occurs
- Each area stores: start_pos, end_pos, center_pos (highest SPC), total_spc, tau

## 5. DTR Deduplication

Before classification, peak areas are deduplicated:
- Areas in the second DTR region are translated to canonical (first region) positions
- Duplicate areas at the same canonical position (within 20bp tolerance) are removed
- Result: unique peak count for classification

**Important:** Canonical positions are used for distance calculation to ensure correct classification (e.g., T7 phage: start at pos 1, end at pos 160 → DTR distance = 159bp).

## 6. Classification

Based on unique peak count after deduplication:

| Start Areas | End Areas | Classification |
|-------------|-----------|----------------|
| 0 | 0 | No_packaging |
| 1 | 0 | PAC |
| 0 | 1 | PAC |
| 1 | 1 | See distance-based rules below |
| >1 or >1 | | Unknown_packaging |

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
- >10% genome → **DTR_outlier_5'/3'**

The suffix (_5' or _3') indicates overhang orientation based on whether end comes before or after start in genomic coordinates.

## 7. Output

- **Mechanism**: Classification string (e.g., "DTR_short_5'")
- **Left termini**: All original start area positions (including both DTR regions)
- **Right termini**: All original end area positions
- **Duplication status**: `true` if all in DTR, `false` if all in ITR, `null` otherwise

## 8. Diagnostic CSV Output

A diagnostic CSV file (`{db_name}_phagetermini.csv`) is generated with one row per position passing the local filter. Columns:

| Column | Description |
|--------|-------------|
| Contig_name | Contig name |
| Sample_name | Sample name |
| Unclipped_ratio | Global unclipped read ratio |
| Coverage_reduced_mean | Global mean coverage_reduced |
| Type | "start" or "end" |
| Position | 1-indexed position |
| SPC | read_starts or read_ends count |
| Clippings | Aggregated clippings within ±max_distance_peaks |
| Observed_ratio | SPC/(SPC+effective_clippings) |
| Primary_reads | Primary reads at position |
| Coverage_reduced | Coverage_reduced at position |
| Tau | SPC/coverage_reduced |
| Kept | "yes" or "no" |
| Discarded_because | Reason if discarded, empty if kept |

Discard reasons:
- `statistical_test: excess_clippings z=X.XX > z_critical=X.XX` - significantly more clippings than expected
