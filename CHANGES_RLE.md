# Adaptive RLE Compression Implementation

## Summary

Replaced the previous outlier-based compression with adaptive run-length encoding (RLE). This simplifies the compression algorithm from 3 parameters to 1, while automatically preserving rapid changes and efficiently compressing flat regions.

## Key Changes

### 1. Database Schema Changes

**Feature Tables** (e.g., `Feature_coverage`, `Feature_reads_starts`, etc.):
- **OLD**: `(Feature_id, Contig_id, Sample_id, Position, Value)`
- **NEW**: `(Feature_id, Contig_id, Sample_id, First_position, Last_position, Value)`

**TempFeatures Table**:
- **OLD**: `(Variable_name, Contig_name, Position, Value)`
- **NEW**: `(Variable_name, Contig_name, First_position, Last_position, Value)`

**Migration Impact**: Existing databases will need to be regenerated. Old databases cannot be read with the new schema.

### 2. Compression Algorithm

**Previous Approach** (3 parameters):
```rust
compress_signal(values, plot_type, step, z_thresh, deriv_thresh, max_points, python_compat)
// - step: regular subsampling interval
// - z_thresh: z-score for value outliers
// - deriv_thresh: z-score for derivative outliers  
// - max_points: hard limit on output
```

**New Approach** (1 parameter):
```rust
compress_signal(values, plot_type, compress_ratio) -> Vec<Run>
// - compress_ratio: relative tolerance (e.g., 0.1 = 10% change threshold)
```

**Algorithm**:
1. Start first run with position 0 and its value
2. For each subsequent position:
   - If `|value - current_run_value| <= compress_ratio * |current_run_value|`, extend run
   - Otherwise, close current run and start new one
3. Track running average of values in each run

**Benefits**:
- **Automatic spike detection**: Rapid changes exceed the threshold and create short runs
- **Efficient flat region compression**: Long constant-coverage stretches become one row
- **Adaptive scaling**: Relative threshold works equally well for high/low coverage regions
- **Simpler**: One intuitive parameter vs three statistical thresholds

### 3. Code Changes

#### Modified Files:

**src/compress.rs**:
- Complete rewrite of compression logic
- New `Run` struct: `{ start_pos, end_pos, value }`
- Added unit tests for RLE behavior

**src/types.rs**:
- `FeaturePoint` struct updated:
  - Removed: `position: i32`
  - Added: `start_pos: i32, end_pos: i32`

**src/processing.rs**:
- `ProcessConfig` simplified:
  - Removed: `step, z_thresh, deriv_thresh, max_points, python_compat`
  - Added: `compress_ratio: f64`
- `add_compressed_feature()` updated to convert `Vec<Run>` to `Vec<FeaturePoint>`

**src/db.rs**:
- Updated table creation DDL for new schema
- Updated `write_features_to_temp_db()` to write runs
- Updated `merge_features()` to read/write runs

**src/main.rs**:
- CLI arguments updated:
  - Removed: `--step, --outlier-threshold, --derivative-threshold, --max-points, --python-compat`
  - Added: `--compress_ratio` (default: 0.1)

**src/lib.rs** (Python bindings):
- Python function signature updated to use `compress_ratio` and `create_indexes`

### 4. Example Usage

**Rust CLI**:
```bash
# Old
mgfeatureviewer -t 4 -g ref.gbk -b bams/ -m coverage,phagetermini -o output.db \
  --step 50 --outlier-threshold 3.0 --derivative-threshold 3.0 --max-points 10000

# New
mgfeatureviewer -t 4 -g ref.gbk -b bams/ -m coverage,phagetermini -o output.db \
  --compress_ratio 0.1
```

**Python bindings**:
```python
# Old
result = mgfeatureviewer_rs.process_all_samples(
    "ref.gbk", "bams/", "output.db", ["coverage"], 4,
    step=50, z_thresh=3.0, deriv_thresh=3.0, max_points=10000
)

# New
result = mgfeatureviewer_rs.process_all_samples(
    "ref.gbk", "bams/", "output.db", ["coverage"], 4,
    compress_ratio=0.1, create_indexes=True
)
```

### 5. Tuning compress_ratio

- **compress_ratio = 0.05** (5%): More sensitive, more runs, captures smaller changes
- **compress_ratio = 0.10** (10%): Balanced (recommended default)
- **compress_ratio = 0.20** (20%): More aggressive compression, only large changes preserved

For most genomic coverage signals, **compress_ratio=0.1** provides good compression while preserving important features.

### 6. Performance Expectations

**Storage Savings**:
- Constant regions: ~1000x compression (1 million positions → ~1000 runs)
- Variable regions: Minimal overhead (similar to original data)
- Overall: Typical 10-100x reduction in database size

**Speed**:
- Compression is O(n) single-pass (faster than old multi-pass outlier detection)
- Database writes reduced proportionally to compression ratio
- Query performance improved (fewer rows to scan)

### 7. Testing

Run unit tests:
```bash
cargo test compress
```

Tests verify:
- Simple RLE behavior with example data
- Constant signals compress to single run
- Spikes create singleton runs
- Empty input handling

## Breaking Changes

⚠️ **This is a breaking change**:
1. Old `.db` files cannot be read (schema incompatible)
2. Python/CLI APIs changed (parameter names/counts)
3. Must regenerate all databases with new version

## Migration Path

To migrate existing data:
1. Keep old BAM files
2. Update code to new version
3. Rerun processing with `--compress_ratio 0.5`
4. Compare outputs to validate (expect ~same visual results, much smaller DB)
