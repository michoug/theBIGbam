use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Float32Array, Int32Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use clap::Parser;
use gb_io::reader::SeqReader;
use indicatif::{ProgressBar, ProgressStyle};
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::htslib;
use rusqlite::Connection;
use std::collections::HashMap;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::sync::Arc;

#[derive(Parser, Debug)]
#[command(name = "mgfeatureviewer_calc")]
#[command(about = "Calculate features from BAM files (Rust implementation)")]
struct Args {
    /// Number of threads
    #[arg(short = 't', long)]
    threads: usize,

    /// Path to GenBank file
    #[arg(short = 'g', long)]
    genbank: PathBuf,

    /// Path to BAM file or directory
    #[arg(short = 'b', long)]
    bam_files: PathBuf,

    /// Modules to compute (comma-separated: coverage,phagetermini,assemblycheck)
    #[arg(short = 'm', long)]
    modules: String,

    /// Output directory
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Annotation tool (optional)
    #[arg(short = 'a', long, default_value = "")]
    annotation_tool: String,

    /// Minimum coverage percentage for contig inclusion
    #[arg(long, default_value = "50")]
    min_coverage: f64,

    /// Step size for compression
    #[arg(long, default_value = "50")]
    step: usize,

    /// Z-score threshold for outliers
    #[arg(long, default_value = "3.0")]
    outlier_threshold: f64,

    /// Derivative threshold for outliers
    #[arg(long, default_value = "3.0")]
    derivative_threshold: f64,

    /// Maximum points after compression
    #[arg(long, default_value = "10000")]
    max_points: usize,
}

#[derive(Clone, Debug)]
struct ContigInfo {
    name: String,
    length: usize,
    annotation_tool: String,
}

#[derive(Clone, Debug)]
struct FeatureAnnotation {
    contig_id: i64,
    start: i64,
    end: i64,
    strand: i64,
    feature_type: String,
    product: Option<String>,
    function: Option<String>,
    phrog: Option<String>,
}

#[derive(Clone, Copy, PartialEq)]
enum SequencingType {
    Long,
    ShortPaired,
    ShortSingle,
}

#[derive(Clone, Copy, PartialEq)]
enum PlotType {
    Curve,
    Bars,
}

struct FeatureConfig {
    plot_type: PlotType,
}

fn get_feature_config(feature: &str) -> FeatureConfig {
    let plot_type = match feature {
        "coverage" | "coverage_reduced" | "read_lengths" | "insert_sizes" => PlotType::Curve,
        _ => PlotType::Bars,
    };
    FeatureConfig { plot_type }
}

/// Detect sequencing type from first N reads
fn detect_sequencing_type(bam_path: &Path) -> Result<SequencingType> {
    let mut bam = bam::Reader::from_path(bam_path)?;
    let mut n_checked = 0;

    for result in bam.records() {
        let record = result?;
        if record.is_unmapped() {
            continue;
        }

        if record.seq_len() > 1000 {
            return Ok(SequencingType::Long);
        }

        if record.is_paired() {
            return Ok(SequencingType::ShortPaired);
        }

        n_checked += 1;
        if n_checked >= 100 {
            break;
        }
    }

    Ok(SequencingType::ShortSingle)
}

/// Read preprocessed data for a single read
struct ReadData {
    ref_start: i64,
    ref_end: i64,
    query_length: i32,
    template_length: i32,
    is_read1: bool,
    is_proper_pair: bool,
    is_reverse: bool,
    cigar: Vec<(u32, u32)>, // (op, len)
    md_tag: Option<Vec<u8>>,
}

/// Process reads for a contig
fn process_reads_for_contig(
    bam: &mut bam::IndexedReader,
    contig_name: &str,
    contig_length: usize,
    modules: &[String],
    _seq_type: SequencingType,
) -> Result<Vec<ReadData>> {
    let tid = bam.header().tid(contig_name.as_bytes());
    if tid.is_none() {
        return Ok(Vec::new());
    }

    bam.fetch((contig_name, 0, contig_length as i64 * 2))?;

    let need_md = modules.contains(&"phagetermini".to_string())
        || modules.contains(&"assemblycheck".to_string());

    let mut reads = Vec::new();

    for result in bam.records() {
        let record = result?;
        if record.is_unmapped() {
            continue;
        }

        let cigar: Vec<(u32, u32)> = record
            .cigar()
            .iter()
            .map(|c| (c.char() as u32, c.len()))
            .collect();

        let md_tag = if need_md {
            record.aux(b"MD").ok().and_then(|aux| {
                match aux {
                    rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes().to_vec()),
                    _ => None,
                }
            })
        } else {
            None
        };

        reads.push(ReadData {
            ref_start: record.pos(),
            ref_end: record.cigar().end_pos(),
            query_length: record.seq_len() as i32,
            template_length: record.insert_size().abs() as i32,
            is_read1: record.is_first_in_template(),
            is_proper_pair: record.is_proper_pair(),
            is_reverse: record.is_reverse(),
            cigar,
            md_tag,
        });
    }

    Ok(reads)
}

/// Calculate coverage feature
fn calculate_coverage(reads: &[ReadData], ref_length: usize) -> Vec<u64> {
    let mut coverage = vec![0u64; ref_length];

    for read in reads {
        let start = (read.ref_start as usize) % ref_length;
        let end = (read.ref_end as usize) % ref_length;

        if start < end {
            for i in start..end {
                coverage[i] += 1;
            }
        } else {
            // Wrap-around for circular contigs
            for i in start..ref_length {
                coverage[i] += 1;
            }
            for i in 0..end {
                coverage[i] += 1;
            }
        }
    }

    coverage
}

/// Check if read starts with a match
fn starts_with_match(cigar: &[(u32, u32)], md: &Option<Vec<u8>>, at_start: bool) -> bool {
    if cigar.is_empty() {
        return false;
    }

    let (op, _) = if at_start { cigar[0] } else { cigar[cigar.len() - 1] };

    // Check for clipping (S=83, H=72) or insertion (I=73)
    let op_char = op as u8 as char;
    if op_char == 'S' || op_char == 'H' || op_char == 'I' {
        return false;
    }

    // Check MD tag
    if let Some(md_bytes) = md {
        if md_bytes.is_empty() {
            return false;
        }
        let val = if at_start { md_bytes[0] } else { md_bytes[md_bytes.len() - 1] };
        // MD starts with number means match
        return val.is_ascii_digit() && val != b'0';
    }

    false
}

/// Calculate phagetermini features
fn calculate_phagetermini(
    reads: &[ReadData],
    ref_length: usize,
    seq_type: SequencingType,
) -> HashMap<String, Vec<u64>> {
    let mut coverage_reduced = vec![0u64; ref_length];
    let mut start_plus = vec![0u64; ref_length];
    let mut start_minus = vec![0u64; ref_length];
    let mut end_plus = vec![0u64; ref_length];
    let mut end_minus = vec![0u64; ref_length];

    for read in reads {
        let check_end = seq_type == SequencingType::Long;

        if starts_with_match(&read.cigar, &read.md_tag, true)
            && (!check_end || starts_with_match(&read.cigar, &read.md_tag, false))
        {
            let start = (read.ref_start as usize) % ref_length;
            let end = (read.ref_end as usize) % ref_length;

            // Update coverage
            if start <= end {
                for i in start..=end.min(ref_length - 1) {
                    coverage_reduced[i] += 1;
                }
            } else {
                for i in start..ref_length {
                    coverage_reduced[i] += 1;
                }
                for i in 0..=end {
                    coverage_reduced[i] += 1;
                }
            }

            // Update starts/ends
            if read.is_reverse {
                start_minus[start] += 1;
                end_minus[end] += 1;
            } else {
                start_plus[start] += 1;
                end_plus[end] += 1;
            }
        }
    }

    // Compute final starts, ends, tau
    let (reads_starts, reads_ends) = match seq_type {
        SequencingType::ShortPaired | SequencingType::ShortSingle => {
            (start_plus, end_minus)
        }
        SequencingType::Long => {
            let starts: Vec<u64> = start_plus.iter().zip(&start_minus).map(|(a, b)| a + b).collect();
            let ends: Vec<u64> = end_plus.iter().zip(&end_minus).map(|(a, b)| a + b).collect();
            (starts, ends)
        }
    };

    let mut results = HashMap::new();
    results.insert("coverage_reduced".to_string(), coverage_reduced);
    results.insert("reads_starts".to_string(), reads_starts);
    results.insert("reads_ends".to_string(), reads_ends);
    // tau is computed later as float

    results
}

/// Calculate assemblycheck features
fn calculate_assemblycheck(
    reads: &[ReadData],
    ref_length: usize,
    seq_type: SequencingType,
) -> HashMap<String, Vec<u64>> {
    let mut results = HashMap::new();

    // Initialize arrays
    let mut left_clippings = vec![0u64; ref_length];
    let mut right_clippings = vec![0u64; ref_length];
    let mut insertions = vec![0u64; ref_length];
    let mut deletions = vec![0u64; ref_length];
    let mut mismatches = vec![0u64; ref_length];

    let mut sum_read_lengths = vec![0u64; ref_length];
    let mut count_read_lengths = vec![0u64; ref_length];
    let mut sum_insert_sizes = vec![0u64; ref_length];
    let mut count_insert_sizes = vec![0u64; ref_length];
    let mut bad_orientations = vec![0u64; ref_length];

    for read in reads {
        let start = (read.ref_start as usize) % ref_length;
        let end = (read.ref_end as usize) % ref_length;

        // Long reads: track read lengths
        if seq_type == SequencingType::Long {
            let positions: Vec<usize> = if start < end {
                (start..end).collect()
            } else {
                (start..ref_length).chain(0..end).collect()
            };
            for pos in positions {
                sum_read_lengths[pos] += read.query_length as u64;
                count_read_lengths[pos] += 1;
            }
        }

        // Short-paired: track insert sizes and bad orientations
        if seq_type == SequencingType::ShortPaired {
            let positions: Vec<usize> = if start < end {
                (start..end).collect()
            } else {
                (start..ref_length).chain(0..end).collect()
            };
            for pos in positions {
                if read.is_read1 && read.template_length > 0 {
                    sum_insert_sizes[pos] += read.template_length as u64;
                    count_insert_sizes[pos] += 1;
                }
                if !read.is_proper_pair {
                    bad_orientations[pos] += 1;
                }
            }
        }

        // Clippings
        if !read.cigar.is_empty() {
            let (first_op, _) = read.cigar[0];
            let (last_op, _) = read.cigar[read.cigar.len() - 1];
            let first_char = first_op as u8 as char;
            let last_char = last_op as u8 as char;

            if first_char == 'S' || first_char == 'H' {
                left_clippings[start] += 1;
            }
            if last_char == 'S' || last_char == 'H' {
                right_clippings[if end > 0 { end - 1 } else { ref_length - 1 }] += 1;
            }
        }

        // Indels from CIGAR
        let mut ref_pos = read.ref_start as usize;
        for (op, len) in &read.cigar {
            let op_char = *op as u8 as char;
            match op_char {
                'I' => {
                    insertions[ref_pos % ref_length] += 1;
                }
                'D' => {
                    for j in 0..(*len as usize) {
                        deletions[(ref_pos + j) % ref_length] += 1;
                    }
                    ref_pos += *len as usize;
                }
                'M' | '=' | 'X' | 'N' => {
                    ref_pos += *len as usize;
                }
                _ => {}
            }
        }

        // Mismatches from MD tag
        if let Some(ref md_bytes) = read.md_tag {
            let mut ref_pos = read.ref_start as usize;
            let mut i = 0;
            while i < md_bytes.len() {
                let c = md_bytes[i];
                if c.is_ascii_digit() {
                    // Parse number
                    let mut num = 0usize;
                    while i < md_bytes.len() && md_bytes[i].is_ascii_digit() {
                        num = num * 10 + (md_bytes[i] - b'0') as usize;
                        i += 1;
                    }
                    ref_pos += num;
                } else if c == b'^' {
                    // Deletion - skip
                    i += 1;
                    while i < md_bytes.len() && md_bytes[i].is_ascii_uppercase() {
                        ref_pos += 1;
                        i += 1;
                    }
                } else if c.is_ascii_uppercase() {
                    // Mismatch
                    mismatches[ref_pos % ref_length] += 1;
                    ref_pos += 1;
                    i += 1;
                } else {
                    i += 1;
                }
            }
        }
    }

    results.insert("left_clippings".to_string(), left_clippings);
    results.insert("right_clippings".to_string(), right_clippings);
    results.insert("insertions".to_string(), insertions);
    results.insert("deletions".to_string(), deletions);
    results.insert("mismatches".to_string(), mismatches);

    if seq_type == SequencingType::Long {
        results.insert("sum_read_lengths".to_string(), sum_read_lengths);
        results.insert("count_read_lengths".to_string(), count_read_lengths);
    }

    if seq_type == SequencingType::ShortPaired {
        results.insert("sum_insert_sizes".to_string(), sum_insert_sizes);
        results.insert("count_insert_sizes".to_string(), count_insert_sizes);
        results.insert("bad_orientations".to_string(), bad_orientations);
    }

    results
}

/// Compress signal for storage
fn compress_signal(
    values: &[f64],
    plot_type: PlotType,
    step: usize,
    z_thresh: f64,
    deriv_thresh: f64,
    max_points: usize,
) -> (Vec<i32>, Vec<f32>) {
    let n = values.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    // Calculate mean and std
    let mean: f64 = values.iter().sum::<f64>() / n as f64;
    let variance: f64 = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    let std = variance.sqrt().max(1e-9);

    // Find value outliers
    let val_outliers: Vec<usize> = values
        .iter()
        .enumerate()
        .filter(|(_, &v)| (v - mean).abs() > z_thresh * std)
        .map(|(i, _)| i)
        .collect();

    let mut keep_idx: Vec<usize> = match plot_type {
        PlotType::Curve => {
            // Regular subsampling
            let mut regular: Vec<usize> = (0..n).step_by(step).collect();

            // Derivative outliers
            let mut derivatives = vec![0.0f64; n];
            for i in 1..n {
                derivatives[i] = values[i] - values[i - 1];
            }
            let deriv_std = {
                let deriv_mean: f64 = derivatives.iter().sum::<f64>() / n as f64;
                let deriv_var: f64 = derivatives.iter().map(|x| (x - deriv_mean).powi(2)).sum::<f64>() / n as f64;
                deriv_var.sqrt().max(1e-9)
            };

            let deriv_outliers: Vec<usize> = derivatives
                .iter()
                .enumerate()
                .filter(|(_, &d)| d.abs() > deriv_thresh * deriv_std)
                .flat_map(|(i, _)| if i > 0 { vec![i - 1, i] } else { vec![i] })
                .filter(|&i| i < n)
                .collect();

            // Combine all indices
            regular.extend(val_outliers);
            regular.extend(deriv_outliers);
            if n > 0 {
                regular.push(n - 1);
            }
            regular.sort_unstable();
            regular.dedup();
            regular
        }
        PlotType::Bars => {
            // Only keep outliers for bars
            val_outliers
        }
    };

    // Apply max_points limit
    if keep_idx.len() > max_points {
        let step_lim = keep_idx.len() / max_points;
        keep_idx = keep_idx.into_iter().step_by(step_lim).collect();
        if let Some(&last) = keep_idx.last() {
            if last != n - 1 && n > 0 {
                keep_idx.push(n - 1);
            }
        }
    }

    // Build output (1-indexed positions)
    let xs: Vec<i32> = keep_idx.iter().map(|&i| (i + 1) as i32).collect();
    let ys: Vec<f32> = keep_idx.iter().map(|&i| values[i] as f32).collect();

    (xs, ys)
}

/// Feature data point
struct FeaturePoint {
    contig_name: String,
    feature: String,
    position: i32,
    value: f32,
}

/// Presence data
struct PresenceData {
    contig_name: String,
    coverage_pct: f32,
}

/// Process a single sample (BAM file)
fn process_sample(
    bam_path: &Path,
    contigs: &[ContigInfo],
    modules: &[String],
    min_coverage: f64,
    step: usize,
    z_thresh: f64,
    deriv_thresh: f64,
    max_points: usize,
) -> Result<(Vec<FeaturePoint>, Vec<PresenceData>, String)> {
    let sample_name = bam_path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .replace("_with_MD", "")
        .to_string();

    let seq_type = detect_sequencing_type(bam_path)?;
    let mut bam = bam::IndexedReader::from_path(bam_path)?;
    // Enable multi-threaded BAM decompression (2 threads per file)
    bam.set_threads(2)?;

    let mut all_features = Vec::new();
    let mut all_presences = Vec::new();

    for contig in contigs {
        // Contig length is halved (doubled for circular mapping)
        let ref_length = contig.length / 2;

        let reads = process_reads_for_contig(&mut bam, &contig.name, ref_length, modules, seq_type)?;

        if reads.is_empty() {
            continue;
        }

        // Calculate coverage percentage
        let mut covered = vec![false; ref_length];
        for read in &reads {
            let start = (read.ref_start.max(0) as usize).min(ref_length);
            let end = (read.ref_end as usize).min(ref_length);
            for i in start..end {
                covered[i] = true;
            }
        }
        let covered_bp: usize = covered.iter().filter(|&&x| x).count();
        let coverage_pct = (covered_bp as f64 / ref_length as f64) * 100.0;

        if coverage_pct < min_coverage {
            continue;
        }

        all_presences.push(PresenceData {
            contig_name: contig.name.clone(),
            coverage_pct: coverage_pct as f32,
        });

        // Calculate features
        if modules.contains(&"coverage".to_string()) || modules.contains(&"assemblycheck".to_string()) {
            let coverage = calculate_coverage(&reads, ref_length);
            let values: Vec<f64> = coverage.iter().map(|&x| x as f64).collect();
            let config = get_feature_config("coverage");
            let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
            for (x, y) in xs.into_iter().zip(ys) {
                all_features.push(FeaturePoint {
                    contig_name: contig.name.clone(),
                    feature: "coverage".to_string(),
                    position: x,
                    value: y,
                });
            }
        }

        if modules.contains(&"phagetermini".to_string()) {
            let pt_features = calculate_phagetermini(&reads, ref_length, seq_type);

            for feature_name in ["coverage_reduced", "reads_starts", "reads_ends"] {
                if let Some(data) = pt_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    let config = get_feature_config(feature_name);
                    let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
                    for (x, y) in xs.into_iter().zip(ys) {
                        all_features.push(FeaturePoint {
                            contig_name: contig.name.clone(),
                            feature: feature_name.to_string(),
                            position: x,
                            value: y,
                        });
                    }
                }
            }

            // Calculate tau
            if let (Some(cov_red), Some(starts), Some(ends)) = (
                pt_features.get("coverage_reduced"),
                pt_features.get("reads_starts"),
                pt_features.get("reads_ends"),
            ) {
                let tau: Vec<f64> = cov_red
                    .iter()
                    .zip(starts)
                    .zip(ends)
                    .map(|((&c, &s), &e)| {
                        if c > 0 {
                            (s + e) as f64 / c as f64
                        } else {
                            0.0
                        }
                    })
                    .collect();
                let config = get_feature_config("tau");
                let (xs, ys) = compress_signal(&tau, config.plot_type, step, z_thresh, deriv_thresh, max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push(FeaturePoint {
                        contig_name: contig.name.clone(),
                        feature: "tau".to_string(),
                        position: x,
                        value: y,
                    });
                }
            }
        }

        if modules.contains(&"assemblycheck".to_string()) {
            let ac_features = calculate_assemblycheck(&reads, ref_length, seq_type);

            // Process each assemblycheck feature
            for feature_name in ["left_clippings", "right_clippings", "insertions", "deletions", "mismatches", "bad_orientations"] {
                if let Some(data) = ac_features.get(feature_name) {
                    let values: Vec<f64> = data.iter().map(|&x| x as f64).collect();
                    let config = get_feature_config(feature_name);
                    let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
                    for (x, y) in xs.into_iter().zip(ys) {
                        all_features.push(FeaturePoint {
                            contig_name: contig.name.clone(),
                            feature: feature_name.to_string(),
                            position: x,
                            value: y,
                        });
                    }
                }
            }

            // Read lengths (long reads)
            if let (Some(sum_rl), Some(count_rl)) = (
                ac_features.get("sum_read_lengths"),
                ac_features.get("count_read_lengths"),
            ) {
                let values: Vec<f64> = sum_rl
                    .iter()
                    .zip(count_rl)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                let config = get_feature_config("read_lengths");
                let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push(FeaturePoint {
                        contig_name: contig.name.clone(),
                        feature: "read_lengths".to_string(),
                        position: x,
                        value: y,
                    });
                }
            }

            // Insert sizes (short-paired)
            if let (Some(sum_is), Some(count_is)) = (
                ac_features.get("sum_insert_sizes"),
                ac_features.get("count_insert_sizes"),
            ) {
                let values: Vec<f64> = sum_is
                    .iter()
                    .zip(count_is)
                    .map(|(&s, &c)| if c > 0 { s as f64 / c as f64 } else { 0.0 })
                    .collect();
                let config = get_feature_config("insert_sizes");
                let (xs, ys) = compress_signal(&values, config.plot_type, step, z_thresh, deriv_thresh, max_points);
                for (x, y) in xs.into_iter().zip(ys) {
                    all_features.push(FeaturePoint {
                        contig_name: contig.name.clone(),
                        feature: "insert_sizes".to_string(),
                        position: x,
                        value: y,
                    });
                }
            }
        }
    }

    Ok((all_features, all_presences, sample_name))
}

/// Write features to Parquet
fn write_features_parquet(features: &[FeaturePoint], path: &Path) -> Result<()> {
    let schema = Schema::new(vec![
        Field::new("contig_name", DataType::Utf8, false),
        Field::new("feature", DataType::Utf8, false),
        Field::new("position", DataType::Int32, false),
        Field::new("value", DataType::Float32, false),
    ]);

    let contig_names: Vec<&str> = features.iter().map(|f| f.contig_name.as_str()).collect();
    let feature_names: Vec<&str> = features.iter().map(|f| f.feature.as_str()).collect();
    let positions: Vec<i32> = features.iter().map(|f| f.position).collect();
    let values: Vec<f32> = features.iter().map(|f| f.value).collect();

    let batch = RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(StringArray::from(contig_names)) as ArrayRef,
            Arc::new(StringArray::from(feature_names)) as ArrayRef,
            Arc::new(Int32Array::from(positions)) as ArrayRef,
            Arc::new(Float32Array::from(values)) as ArrayRef,
        ],
    )?;

    let file = File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}

/// Write presences to Parquet
fn write_presences_parquet(presences: &[PresenceData], path: &Path) -> Result<()> {
    let schema = Schema::new(vec![
        Field::new("contig_name", DataType::Utf8, false),
        Field::new("coverage_pct", DataType::Float32, false),
    ]);

    let contig_names: Vec<&str> = presences.iter().map(|p| p.contig_name.as_str()).collect();
    let coverage_pcts: Vec<f32> = presences.iter().map(|p| p.coverage_pct).collect();

    let batch = RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(StringArray::from(contig_names)) as ArrayRef,
            Arc::new(Float32Array::from(coverage_pcts)) as ArrayRef,
        ],
    )?;

    let file = File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}

/// Create metadata SQLite database
fn create_metadata_db(
    db_path: &Path,
    contigs: &[ContigInfo],
    annotations: &[FeatureAnnotation],
) -> Result<()> {
    let conn = Connection::open(db_path)?;

    // Create tables
    conn.execute(
        "CREATE TABLE Contig (
            Contig_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_name TEXT UNIQUE,
            Contig_length INTEGER,
            Annotation_tool TEXT
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Sample (
            Sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Sample_name TEXT UNIQUE
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Presences (
            Presence_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER,
            Sample_id INTEGER,
            Coverage_percentage REAL,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Contig_annotation (
            Contig_annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER,
            Start INTEGER,
            End INTEGER,
            Strand INTEGER,
            Type TEXT,
            Product TEXT,
            Function TEXT,
            Phrog INTEGER,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id)
        )",
        [],
    )?;

    conn.execute(
        "CREATE TABLE Variable (
            Variable_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Variable_name TEXT UNIQUE,
            Subplot TEXT,
            Module TEXT,
            Type TEXT,
            Color TEXT,
            Alpha REAL,
            Fill_alpha REAL,
            Size REAL,
            Title TEXT,
            Help TEXT,
            Feature_table_name TEXT
        )",
        [],
    )?;

    // Insert contigs
    for contig in contigs {
        conn.execute(
            "INSERT INTO Contig (Contig_name, Contig_length, Annotation_tool) VALUES (?1, ?2, ?3)",
            [&contig.name, &(contig.length / 2).to_string(), &contig.annotation_tool],
        )?;
    }

    // Insert annotations in a transaction for speed
    conn.execute("BEGIN TRANSACTION", [])?;
    for ann in annotations {
        conn.execute(
            "INSERT INTO Contig_annotation (Contig_id, Start, End, Strand, Type, Product, Function, Phrog)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            rusqlite::params![
                ann.contig_id,
                ann.start,
                ann.end,
                ann.strand,
                ann.feature_type,
                ann.product,
                ann.function,
                ann.phrog,
            ],
        )?;
    }
    conn.execute("COMMIT", [])?;

    // Insert default variables (matching Python)
    let variables = [
        ("coverage", "Coverage", "Coverage", "curve", "#333333", 0.8, 0.4, 1.0, "Coverage depth", ""),
        ("coverage_reduced", "Coverage reduced", "Phage termini", "curve", "#00c53b", 0.8, 0.4, 1.0, "Coverage reduced", ""),
        ("reads_starts", "Reads termini", "Phage termini", "bars", "#215732", 0.6, 0.4, 1.0, "Read Starts", ""),
        ("reads_ends", "Reads termini", "Phage termini", "bars", "#6cc24a", 0.6, 0.4, 1.0, "Read Ends", ""),
        ("tau", "Tau", "Phage termini", "bars", "#44883e", 0.6, 0.4, 1.0, "Tau", ""),
        ("read_lengths", "Read lengths", "Assembly check", "curve", "#ed8b00", 0.8, 0.4, 1.0, "Read Lengths", ""),
        ("insert_sizes", "Insert sizes", "Assembly check", "curve", "#ed8b00", 0.8, 0.4, 1.0, "Insert Sizes", ""),
        ("bad_orientations", "Bad orientations", "Assembly check", "bars", "#c94009", 0.6, 0.4, 1.0, "Bad Orientations", ""),
        ("left_clippings", "Clippings", "Assembly check", "bars", "#7f0091", 0.6, 0.4, 1.0, "Left Clippings", ""),
        ("right_clippings", "Clippings", "Assembly check", "bars", "#8e43e7", 0.6, 0.4, 1.0, "Right Clippings", ""),
        ("insertions", "Indels", "Assembly check", "bars", "#e50001", 0.6, 0.4, 1.0, "Insertions", ""),
        ("deletions", "Indels", "Assembly check", "bars", "#97011a", 0.6, 0.4, 1.0, "Deletions", ""),
        ("mismatches", "Mismatches", "Assembly check", "bars", "#5a0f0b", 0.6, 0.4, 1.0, "Mismatches", ""),
    ];

    for (name, subplot, module, vtype, color, alpha, fill_alpha, size, title, help) in variables {
        let table_name = format!("Feature_{}", name);
        conn.execute(
            "INSERT INTO Variable (Variable_name, Subplot, Module, Type, Color, Alpha, Fill_alpha, Size, Title, Help, Feature_table_name)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            rusqlite::params![name, subplot, module, vtype, color, alpha, fill_alpha, size, title, help, table_name],
        )?;
    }

    Ok(())
}

/// Check if location is complemented (reverse strand)
fn is_complement(loc: &gb_io::seq::Location) -> bool {
    matches!(loc, gb_io::seq::Location::Complement(_))
}

/// Parse GenBank file
fn parse_genbank(path: &Path, annotation_tool: &str) -> Result<(Vec<ContigInfo>, Vec<FeatureAnnotation>)> {
    let file = File::open(path).context("Failed to open GenBank file")?;
    let reader = SeqReader::new(file);

    let mut contigs = Vec::new();
    let mut annotations = Vec::new();
    let mut contig_id = 1i64;

    for result in reader {
        let seq = result.context("Failed to parse GenBank record")?;

        let name = seq.name.clone().unwrap_or_else(|| format!("contig_{}", contig_id));
        let length = seq.seq.len();

        contigs.push(ContigInfo {
            name: name.clone(),
            length,
            annotation_tool: annotation_tool.to_string(),
        });

        // Extract features
        for feature in &seq.features {
            let feature_type = feature.kind.to_string();
            if feature_type == "source" || feature_type == "gene" {
                continue;
            }

            // Use find_bounds() for start/end (returns 0-based exclusive range)
            let (start, end) = match feature.location.find_bounds() {
                Ok((s, e)) => (s + 1, e),  // Convert to 1-based
                Err(_) => continue,
            };

            // Check if complement (reverse strand)
            let strand = if is_complement(&feature.location) { -1 } else { 1 };

            let product = feature.qualifier_values("product").next().map(|s| s.to_string());
            let function = feature.qualifier_values("function").next().map(|s| s.to_string());
            let phrog = feature.qualifier_values("phrog").next().map(|s| s.to_string());

            annotations.push(FeatureAnnotation {
                contig_id,
                start,
                end,
                strand,
                feature_type,
                product,
                function,
                phrog,
            });
        }

        contig_id += 1;
    }

    Ok((contigs, annotations))
}

/// Update Sample and Presences tables in SQLite database
fn update_sample_presences(
    db_path: &Path,
    results: &[Result<(String, Vec<PresenceData>, f64)>],
    contigs: &[ContigInfo],
) -> Result<()> {
    let conn = Connection::open(db_path)?;

    // Build contig name -> id mapping
    let contig_name_to_id: HashMap<String, i64> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.clone(), (i + 1) as i64))
        .collect();

    // Insert samples and presences
    conn.execute("BEGIN TRANSACTION", [])?;

    let mut sample_id = 1i64;
    for result in results {
        if let Ok((sample_name, presences, _)) = result {
            // Insert sample
            conn.execute(
                "INSERT INTO Sample (Sample_name) VALUES (?1)",
                [sample_name],
            )?;

            // Insert presences for this sample
            for presence in presences {
                if let Some(&contig_id) = contig_name_to_id.get(&presence.contig_name) {
                    conn.execute(
                        "INSERT INTO Presences (Contig_id, Sample_id, Coverage_percentage) VALUES (?1, ?2, ?3)",
                        rusqlite::params![contig_id, sample_id, presence.coverage_pct],
                    )?;
                }
            }
            sample_id += 1;
        }
    }

    conn.execute("COMMIT", [])?;
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Suppress htslib warnings about index file timestamps
    unsafe {
        htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
    }

    // Set up rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    println!("### Parsing GenBank file...");
    let (contigs, annotations) = parse_genbank(&args.genbank, &args.annotation_tool)?;
    println!("Found {} contigs with {} annotations", contigs.len(), annotations.len());

    // Get BAM files
    let bam_files: Vec<PathBuf> = if args.bam_files.is_dir() {
        walkdir::WalkDir::new(&args.bam_files)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.path().extension().map(|x| x == "bam").unwrap_or(false))
            .map(|e| e.path().to_path_buf())
            .collect()
    } else {
        vec![args.bam_files.clone()]
    };

    if bam_files.is_empty() {
        anyhow::bail!("No BAM files found");
    }

    // Sort by file size (largest first) for load balancing
    let mut bam_files = bam_files;
    bam_files.sort_by_key(|p| std::cmp::Reverse(fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    println!("Found {} BAM files", bam_files.len());

    // Create output directories
    if args.output.exists() {
        anyhow::bail!("Output directory already exists: {:?}", args.output);
    }
    fs::create_dir_all(&args.output)?;
    fs::create_dir_all(args.output.join("features"))?;
    fs::create_dir_all(args.output.join("presences"))?;

    // Sync to ensure directories are created before SQLite
    std::thread::sleep(std::time::Duration::from_millis(100));

    // Create metadata database - use canonical path
    let output_dir = fs::canonicalize(&args.output)?;
    let db_path = output_dir.join("metadata.db");
    create_metadata_db(&db_path, &contigs, &annotations)?;

    // Parse modules
    let modules: Vec<String> = args.modules.split(',').map(|s| s.trim().to_string()).collect();

    println!("\n### Processing {} samples with {} threads", bam_files.len(), args.threads);
    println!("Modules: {}", modules.join(", "));
    println!();

    // Create progress bar
    let pb = ProgressBar::new(bam_files.len() as u64);
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let start_time = std::time::Instant::now();
    use std::sync::atomic::{AtomicUsize, Ordering};
    let error_count = AtomicUsize::new(0);
    let completed_count = AtomicUsize::new(0);
    let total = bam_files.len();

    // Process samples in parallel - collect presences for SQLite
    let results: Vec<Result<(String, Vec<PresenceData>, f64)>> = bam_files
        .par_iter()
        .map(|bam_path| {
            let sample_start = std::time::Instant::now();

            let result = process_sample(
                bam_path,
                &contigs,
                &modules,
                args.min_coverage,
                args.step,
                args.outlier_threshold,
                args.derivative_threshold,
                args.max_points,
            );

            match result {
                Ok((features, presences, sample_name)) => {
                    // Write Parquet files
                    if !features.is_empty() {
                        let features_path = args.output.join("features").join(format!("{}.parquet", sample_name));
                        write_features_parquet(&features, &features_path)?;
                    }

                    if !presences.is_empty() {
                        let presences_path = args.output.join("presences").join(format!("{}.parquet", sample_name));
                        write_presences_parquet(&presences, &presences_path)?;
                    }

                    let sample_time = sample_start.elapsed().as_secs_f64();
                    let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
                    pb.set_message(format!("[{}/{}] {} ({:.2}s)", done, total, sample_name, sample_time));
                    pb.inc(1);
                    Ok((sample_name, presences, sample_time))
                }
                Err(e) => {
                    let sample_time = sample_start.elapsed().as_secs_f64();
                    let done = completed_count.fetch_add(1, Ordering::SeqCst) + 1;
                    pb.set_message(format!("[{}/{}] ERR: {} ({:.2}s)", done, total, bam_path.file_name().unwrap_or_default().to_string_lossy(), sample_time));
                    pb.inc(1);
                    error_count.fetch_add(1, Ordering::Relaxed);
                    Err(e)
                }
            }
        })
        .collect();

    pb.finish_with_message("Done");

    // Insert Sample and Presences into SQLite
    println!("Updating metadata database...");
    update_sample_presences(&db_path, &results, &contigs)?;

    let elapsed = start_time.elapsed();
    let total_samples = bam_files.len();
    let successful = total_samples - error_count.load(Ordering::Relaxed);

    // Calculate statistics
    let sample_times: Vec<f64> = results.iter().filter_map(|r| r.as_ref().ok().map(|(_, _, t)| *t)).collect();
    let avg_time = if !sample_times.is_empty() {
        sample_times.iter().sum::<f64>() / sample_times.len() as f64
    } else {
        0.0
    };

    println!();
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("### Complete");
    println!();
    println!("  Samples processed: {}/{}", successful, total_samples);
    println!("  Total time:        {:.2}s", elapsed.as_secs_f64());
    println!("  Avg time/sample:   {:.2}s", avg_time);
    println!();
    println!("  Output: {:?}", args.output);
    println!("    - metadata.db");
    println!("    - features/*.parquet ({} files)", successful);
    println!("    - presences/*.parquet ({} files)", successful);
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    if error_count.load(Ordering::Relaxed) > 0 {
        eprintln!("\nWarning: {} samples failed to process", error_count.load(Ordering::Relaxed));
    }

    Ok(())
}
