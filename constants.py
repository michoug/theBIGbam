# Define function-to-color mapping
# Use the color scheme from pharokka
PHAROKKA_COLORS = {
    "vfdb_card": "#FF0000",
    "unknown function": "#AAAAAA",
    "other": "#4deeea",
    "tail": "#74ee15",
    "transcription regulation": "#ffe700",
    "dna, rna and nucleotide metabolism": "#f000ff",
    "lysis": "#001eff",
    "moron, auxiliary metabolic gene and host takeover": "#8900ff",
    "integration and excision": "#E0B0FF",
    "head and packaging": "#ff008d",
    "connector": "#5A5A5A",
}

# Define subplots characteristics
def config_feature_subplot(plot_type, color, alpha, size, title):
    return {
        "type_picked": plot_type,
        "color_picked": color,
        "alpha_picked": alpha,
        "size_picked": size,
        "title_picked": title
    }

FEATURE_SUBPLOTS = {
    "coverage": config_feature_subplot("curve", "black", 0.5, 2, "Coverage Depth"),
    "coverage_reduced": config_feature_subplot("curve", "black", 0.5, 2, "Coverage Depth (only reads starting and ending with a match)"),
    "reads_starts": config_feature_subplot("dots", "blue", 0.3, 6, "Reads' Starts (more than deviation_factor*std away from the mean)"),
    "reads_ends": config_feature_subplot("dots", "red", 0.3, 6, "Reads' Ends (more than deviation_factor*std away from the mean)"),
    "tau": config_feature_subplot("dots", "green", 0.3, 6, "Tau (more than deviation_factor*std away from the mean)"),
}