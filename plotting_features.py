import constants
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
import numpy as np
import pysam

# Function dedicated to specific features
def reduce_position(pos, ref_length):
    return pos % ref_length

def calculate_coverage(coverage, read, ref_length):
    for block_start, block_end in read.get_blocks():
        if block_start < ref_length and block_end > ref_length:  # circular boundary
            coverage[block_start:ref_length] += 1
            coverage[0:reduce_position(block_end, ref_length)] += 1
        else:
            coverage[reduce_position(block_start, ref_length):reduce_position(block_end, ref_length)] += 1

    return coverage

def starts_with_match(read, start):
    """
    Return True if the first aligned base is a match (not clipped, not insertion, not mismatch).
    """
    if read.is_unmapped or read.cigartuples is None:
        return False

    # Check clipping at start
    first_op, _ = read.cigartuples[0]
    if not start:
        first_op, _ = read.cigartuples[-1]  # check end if start is False
    if first_op in (4, 5):  # soft or hard clip
        return False
    if first_op == 1:  # insertion
        return False

    # Now check MD tag
    try:
        md = read.get_tag("MD")
    except KeyError:
        raise ValueError("No MD tag present in BAM")

    md_status = md[0].isdigit()
    if not start:
        md_status = md[-1].isdigit()  # check end if start is False
    return md_status

def calculate_reads_starts_and_ends(read, temporary_feature_dict, ref_length):
    if temporary_feature_dict:
        if starts_with_match(read, start=True) and starts_with_match(read, start=False):
            start = read.reference_start
            end = read.reference_end

            if start < ref_length and end > ref_length:
                temporary_feature_dict["coverage_reduced"][start:ref_length] += 1
                temporary_feature_dict["coverage_reduced"][0:reduce_position(end, ref_length)] += 1
            else:
                temporary_feature_dict["coverage_reduced"][reduce_position(start, ref_length):reduce_position(end, ref_length)] += 1

            start_pos = reduce_position(read.reference_start, ref_length)
            end_pos = reduce_position(read.reference_end - 1, ref_length)  # end is exclusive

            if read.is_reverse:
                temporary_feature_dict["start_minus"][start_pos] += 1
                temporary_feature_dict["end_minus"][end_pos] += 1
            else:
                temporary_feature_dict["start_plus"][start_pos] += 1
                temporary_feature_dict["end_plus"][end_pos] += 1

    return temporary_feature_dict

def compute_final_starts_ends_and_tau(feature_dict, temporary_feature_dict, sequencing_type, ref_length):
    # Combine per sequencing_type
    feature_dict["coverage_reduced"] = temporary_feature_dict["coverage_reduced"]
    if sequencing_type == "short":
        feature_dict["reads_starts"] = temporary_feature_dict["start_plus"]
        feature_dict["reads_ends"] = temporary_feature_dict["end_minus"]
    elif sequencing_type == "long":
        feature_dict["reads_starts"] = temporary_feature_dict["start_plus"] + temporary_feature_dict["start_minus"]
        feature_dict["reads_ends"] = temporary_feature_dict["end_plus"] + temporary_feature_dict["end_minus"]
    else:
        raise ValueError("Sequencing_type must be 'short' or 'long'")

    # --- COMPUTE TAU ---
    tau = np.zeros(ref_length, dtype=np.float32)
    mask = feature_dict["coverage_reduced"] > 0
    tau[mask] = (feature_dict["reads_starts"][mask] + feature_dict["reads_ends"][mask]) / feature_dict["coverage_reduced"][mask]
    feature_dict["tau"] = tau

    return feature_dict

### Main logic to get features
def get_features(bamfile, features, reference, ref_length, sequencing_type):    
    # Initialize arrays for all requested features
    feature_dict = {feature: np.zeros(ref_length, dtype=np.uint64) for feature in features}
    
    # Temporary array for reduced coverage (reads starting/ending with a match)
    temporary_starts_ends_dict = {}
    if {"reads_starts", "reads_ends"}.intersection(features) or "tau" in features:
        temporary_feature = ["start_plus", "start_minus", "end_plus", "end_minus", "coverage_reduced"]
        temporary_starts_ends_dict = {feature: np.zeros(ref_length, dtype=np.uint64) for feature in temporary_feature}

    for read in bamfile.fetch(reference):
        if read.is_unmapped:
            continue

        # --- COVERAGE ---
        feature_dict["coverage"] = calculate_coverage(feature_dict["coverage"], read, ref_length) if "coverage" in features else feature_dict["coverage"]
        # --- START/END BY STRAND ---
        temporary_starts_ends_dict = calculate_reads_starts_and_ends(read, temporary_starts_ends_dict, ref_length) if temporary_starts_ends_dict else temporary_starts_ends_dict

    # --- FINALIZE START/END AND TAU ---
    feature_dict = compute_final_starts_ends_and_tau(feature_dict, temporary_starts_ends_dict, sequencing_type, ref_length) if temporary_starts_ends_dict else feature_dict

    return feature_dict

def smooth_values(feature_values, ref_length, window_size):
    # Aggregate by window
    n_windows = (ref_length + window_size - 1) // window_size
    window_cov = np.add.reduceat(feature_values, np.arange(0, ref_length, window_size))
    # If genome length not divisible by window size, trim
    window_cov = window_cov[:n_windows]

    # Convert to dict of window -> coverage_depth (sum or mean)
    coverage_dict = {
        i + 1: float(window_cov[i]) / window_size
        for i in range(n_windows)
    }

    return coverage_dict

### Plotting functions
def make_bokeh_subplot(feature, xx, yy, width, height, x_range):
    type_picked = constants.FEATURE_SUBPLOTS[feature]["type_picked"]
    color_picked = constants.FEATURE_SUBPLOTS[feature]["color_picked"]
    alpha_picked = constants.FEATURE_SUBPLOTS[feature]["alpha_picked"]
    size_picked = constants.FEATURE_SUBPLOTS[feature]["size_picked"]
    title_picked = constants.FEATURE_SUBPLOTS[feature]["title_picked"]

    p = figure(
        width=width,
        height=height,
        title=title_picked,
        x_range=x_range,
        tools="xpan,xwheel_zoom,reset,save"
    )

    # Part specific to the type of subplot
    if type_picked == "curve":
        p.line(
            x=xx,
            y=yy,
            line_color=color_picked,
            line_alpha=alpha_picked,
            line_width=size_picked,
        )
    elif type_picked == "dots":
        p.scatter(
            x=xx,
            y=yy,
            color=color_picked,
            alpha=alpha_picked,
            size=size_picked
        )

    # Add hover
    source = ColumnDataSource(data=dict(x=xx, y=yy))
    p.line(x='x', y='y', source=source, line_alpha=0, line_width=1)  # line_alpha=0 makes it invisible
    hover = HoverTool(tooltips=[("Position", "@x"), (title_picked, "@y")], mode='vline')
    p.add_tools(hover)

    # A clean style like your matplotlib setup
    p.xgrid.visible = False

    p.y_range.start = 0
    p.yaxis.axis_label = title_picked
    p.yaxis.axis_label_text_font_size = "10pt"
    p.yaxis.axis_label_standoff = 0
    p.ygrid.grid_line_alpha = 0.2
    p.yaxis.axis_label = None
    
    p.outline_line_color = None  # hides top/right borders
    p.min_border_left = 40
    p.min_border_right = 10

    return p

def prepare_subplot(feature, feature_values, locus_size, max_visible_width, subplot_size, shared_xrange, window_size, deviation_factor = 3):
    feature_values_averaged = smooth_values(feature_values, locus_size, window_size)

    # Define window midpoints and y-values
    xx = np.arange(window_size / 2, locus_size, window_size)
    yy = np.array([feature_values_averaged[i + 1] for i in range(len(xx))], dtype=float)

    type_picked = constants.FEATURE_SUBPLOTS[feature]["type_picked"]
    if type_picked == "dots":
        # --- Filter values that diverge strongly ---
        mean_y = np.mean(yy)
        std_y = np.std(yy)
        mask = np.abs(yy - mean_y) > (deviation_factor * std_y)
        xx = xx[mask]
        yy = yy[mask]

    feature_subplot = make_bokeh_subplot(feature, xx, yy, max_visible_width, subplot_size, shared_xrange)
    return feature_subplot

### One function to rule them all
def adding_subplots(mapping_file, features_list, locus_name, locus_size, sequencing_type, max_visible_width, subplot_size, shared_xrange, window_size):
    # Read bam once to calculate all features
    print("Read bam once to calculate all features...", flush=True)
    bam_file = pysam.AlignmentFile(mapping_file, "rb")
    feature_values = get_features(bam_file, features_list, locus_name, locus_size, sequencing_type)

    print("Plotting feature subplots...", flush=True)
    subplots = []
    for feature in features_list:
        subplot_feature = prepare_subplot(feature, feature_values[feature], locus_size, max_visible_width, subplot_size, shared_xrange, window_size)
        subplots.append(subplot_feature)

    bam_file.close()
    return subplots