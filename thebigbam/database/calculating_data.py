import argparse, sys, os
import duckdb
from pathlib import Path
from multiprocessing import cpu_count

# Import Rust bindings (required)
try:
    import thebigbam_rs as _rust
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _rust = None


def _store_contig_sequences(db_path, assembly_path=None, genbank_path=None):
    """Store contig sequences in the database for sequence visualization.

    Reads sequences from the assembly FASTA or GenBank file and stores them
    in a Contig_sequence table, matched by contig name.

    Args:
        db_path: Path to the output DuckDB database
        assembly_path: Path to assembly FASTA file (preferred source)
        genbank_path: Path to GenBank annotation file (fallback source; GFF3 has no sequences)
    """
    from Bio import SeqIO

    # Pick the source: assembly FASTA if available, else GenBank
    source_path = None
    source_format = None
    if assembly_path and os.path.exists(assembly_path):
        source_path = assembly_path
        source_format = "fasta"
    elif genbank_path and os.path.exists(genbank_path):
        # Skip GFF3 files — they don't contain sequences
        if genbank_path.lower().endswith(('.gff', '.gff3')):
            print("  Skipping sequence storage: GFF3 files don't contain sequences.", flush=True)
            return
        source_path = genbank_path
        source_format = "genbank"
    else:
        return

    # Read sequences into a dict
    sequences = {}
    try:
        for record in SeqIO.parse(source_path, source_format):
            sequences[record.id] = str(record.seq)
    except Exception as e:
        print(f"  WARNING: Could not read sequences from {source_path}: {e}", flush=True)
        return

    if not sequences:
        print("  No sequences found in source file.", flush=True)
        return

    # Open the DB and store sequences matched to existing contigs
    try:
        db_conn = duckdb.connect(db_path)
        cur = db_conn.cursor()

        cur.execute("""
            CREATE TABLE IF NOT EXISTS Contig_sequence (
                Contig_id INTEGER PRIMARY KEY,
                Sequence TEXT NOT NULL
            )
        """)

        # Get existing contig names and IDs
        cur.execute("SELECT Contig_id, Contig_name FROM Contig")
        contig_rows = cur.fetchall()

        inserted = 0
        for contig_id, contig_name in contig_rows:
            if contig_name in sequences:
                cur.execute(
                    "INSERT INTO Contig_sequence (Contig_id, Sequence) VALUES (?, ?)",
                    (contig_id, sequences[contig_name])
                )
                inserted += 1

        db_conn.close()
        print(f"  Stored sequences for {inserted}/{len(contig_rows)} contigs.", flush=True)
    except Exception as e:
        print(f"  WARNING: Could not store sequences in database: {e}", flush=True)


def _create_codon_table_and_translated_view(db_path):
    """Create Codon_table and Contig_annotation_translated table.

    Translates CDS annotations to protein sequences using the standard genetic code.
    Only runs if both Contig_sequence and Contig_annotation tables exist.

    Args:
        db_path: Path to the DuckDB database
    """
    from Bio.Seq import Seq

    try:
        db_conn = duckdb.connect(db_path)
        cur = db_conn.cursor()

        # Check that both required tables exist
        cur.execute("SELECT table_name FROM information_schema.tables WHERE table_name IN ('Contig_sequence', 'Contig_annotation')")
        existing = {r[0] for r in cur.fetchall()}
        if not {'Contig_sequence', 'Contig_annotation'}.issubset(existing):
            db_conn.close()
            return

        # --- Create and populate Codon_table ---
        cur.execute("""
            CREATE TABLE IF NOT EXISTS Codon_table (
                Codon TEXT PRIMARY KEY,
                AminoAcid TEXT NOT NULL,
                AminoAcid_name TEXT NOT NULL,
                Color TEXT NOT NULL
            )
        """)

        # Standard genetic code: codon -> (1-letter AA, full name, color)
        # Colors by biochemical property (chosen to avoid nucleotide colors:
        #   A=#d62728, T=#2ca02c, G=#ff7f0e, C=#1f77b4):
        #   Hydrophobic aliphatic (G,A,V,L,I,P): #2d6a4f (dark teal-green)
        #   Aromatic (F,W,Y): #b5838d (muted rose)
        #   Polar uncharged (S,T,N,Q,C,M): #457b9d (steel blue)
        #   Positively charged (K,R,H): #9b2226 (dark crimson)
        #   Negatively charged (D,E): #e9c46a (golden yellow)
        #   Stop (*): #6a3d9a (purple)
        aa_info = {
            'G': ('Glycine', '#2d6a4f'), 'A': ('Alanine', '#2d6a4f'),
            'V': ('Valine', '#2d6a4f'), 'L': ('Leucine', '#2d6a4f'),
            'I': ('Isoleucine', '#2d6a4f'), 'P': ('Proline', '#2d6a4f'),
            'F': ('Phenylalanine', '#b5838d'), 'W': ('Tryptophan', '#b5838d'),
            'Y': ('Tyrosine', '#b5838d'),
            'S': ('Serine', '#457b9d'), 'T': ('Threonine', '#457b9d'),
            'N': ('Asparagine', '#457b9d'), 'Q': ('Glutamine', '#457b9d'),
            'C': ('Cysteine', '#457b9d'), 'M': ('Methionine', '#457b9d'),
            'K': ('Lysine', '#9b2226'), 'R': ('Arginine', '#9b2226'),
            'H': ('Histidine', '#9b2226'),
            'D': ('Aspartate', '#e9c46a'), 'E': ('Glutamate', '#e9c46a'),
            '*': ('Stop', '#6a3d9a'),
        }

        # Translate all 64 codons
        bases = 'TCAG'
        codon_rows = []
        for b1 in bases:
            for b2 in bases:
                for b3 in bases:
                    codon = b1 + b2 + b3
                    aa = str(Seq(codon).translate())
                    name, color = aa_info.get(aa, ('Unknown', '#999999'))
                    codon_rows.append((codon, aa, name, color))

        cur.execute("DELETE FROM Codon_table")
        cur.executemany(
            "INSERT INTO Codon_table (Codon, AminoAcid, AminoAcid_name, Color) VALUES (?, ?, ?, ?)",
            codon_rows
        )

        # --- Create and populate Contig_annotation_translated ---
        cur.execute("DROP TABLE IF EXISTS Contig_annotation_translated")

        # Get column names from Contig_annotation
        cur.execute("SELECT column_name FROM information_schema.columns WHERE table_name = 'Contig_annotation' ORDER BY ordinal_position")
        ca_columns = [r[0] for r in cur.fetchall()]
        col_list = ', '.join(f'"{c}"' for c in ca_columns)

        cur.execute(f"""
            CREATE TABLE Contig_annotation_translated AS
            SELECT {col_list},
                CAST(NULL AS TEXT) AS Nucleotide_sequence,
                CAST(NULL AS TEXT) AS Protein_sequence
            FROM Contig_annotation WHERE 1=0
        """)

        # Get all annotations
        cur.execute(f"SELECT {col_list} FROM Contig_annotation")
        annotations = cur.fetchall()

        # Get all contig sequences keyed by Contig_id
        cur.execute("SELECT Contig_id, Sequence FROM Contig_sequence")
        contig_seqs = {r[0]: r[1] for r in cur.fetchall()}

        # Find column indices
        col_idx = {name: i for i, name in enumerate(ca_columns)}
        type_idx = col_idx.get('Type')
        start_idx = col_idx.get('Start')
        end_idx = col_idx.get('End')
        strand_idx = col_idx.get('Strand')
        contig_id_idx = col_idx.get('Contig_id')

        if any(idx is None for idx in [type_idx, start_idx, end_idx, strand_idx, contig_id_idx]):
            print("  WARNING: Contig_annotation missing required columns for translation.", flush=True)
            db_conn.close()
            return

        insert_rows = []
        translated_count = 0
        for row in annotations:
            row_list = list(row)
            feat_type = row_list[type_idx]
            if feat_type == 'CDS':
                contig_id = row_list[contig_id_idx]
                seq_full = contig_seqs.get(contig_id)
                if seq_full:
                    start = int(row_list[start_idx])
                    end = int(row_list[end_idx])
                    strand = int(row_list[strand_idx]) if row_list[strand_idx] is not None else 1
                    # Extract DNA subsequence (1-based inclusive coordinates)
                    dna = seq_full[start - 1:end]
                    if strand == -1:
                        dna = str(Seq(dna).reverse_complement())
                    # Translate (allow partial codons at end)
                    try:
                        protein = str(Seq(dna).translate())
                    except Exception:
                        protein = None
                        dna = None
                    row_list.append(dna)
                    row_list.append(protein)
                    if protein is not None:
                        translated_count += 1
                else:
                    row_list.append(None)
                    row_list.append(None)
            else:
                row_list.append(None)
                row_list.append(None)
            insert_rows.append(tuple(row_list))

        # Build placeholders for insert
        placeholders = ', '.join(['?'] * (len(ca_columns) + 2))
        all_cols = col_list + ', "Nucleotide_sequence", "Protein_sequence"'
        cur.executemany(
            f"INSERT INTO Contig_annotation_translated ({all_cols}) VALUES ({placeholders})",
            insert_rows
        )

        db_conn.close()
        print(f"  Translated {translated_count} CDS annotations to protein sequences.", flush=True)

    except Exception as e:
        print(f"  WARNING: Could not create translated annotations: {e}", flush=True)


def _create_mismatches_with_codons(db_path):
    """Create Feature_mismatches_with_codons table with synonymous/non-synonymous annotation.

    For each mismatch in a CDS, computes the codon change and whether the amino acid
    changes (non-synonymous) or stays the same (synonymous). Positions not in any CDS
    remain NULL (intergenic).

    Args:
        db_path: Path to the DuckDB database
    """
    import bisect

    try:
        db_conn = duckdb.connect(db_path)
        cur = db_conn.cursor()

        # Check that both required tables exist
        cur.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_name IN ('Feature_mismatches', 'Contig_annotation_translated')"
        )
        existing = {r[0] for r in cur.fetchall()}
        if not {'Feature_mismatches', 'Contig_annotation_translated'}.issubset(existing):
            db_conn.close()
            return

        # Create the enriched table as a copy of Feature_mismatches + 3 new columns
        cur.execute("DROP TABLE IF EXISTS Feature_mismatches_with_codons")
        cur.execute("""
            CREATE TABLE Feature_mismatches_with_codons AS
            SELECT *, CAST(NULL AS TEXT) AS Codon_category,
                      CAST(NULL AS TEXT) AS Codon_change,
                      CAST(NULL AS TEXT) AS AA_change
            FROM Feature_mismatches
        """)

        # Build codon -> amino acid lookup from Codon_table
        cur.execute("SELECT Codon, AminoAcid, AminoAcid_name FROM Codon_table")
        codon_to_aa = {}       # codon -> 1-letter AA
        codon_to_name = {}     # codon -> full name
        for codon, aa, aa_name in cur.fetchall():
            codon_to_aa[codon] = aa
            codon_to_name[codon] = aa_name

        if not codon_to_aa:
            db_conn.close()
            return

        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        # Get all distinct contig IDs that have mismatches
        cur.execute("SELECT DISTINCT Contig_id FROM Feature_mismatches_with_codons WHERE Sequence IS NOT NULL AND Sequence != ''")
        contig_ids = [r[0] for r in cur.fetchall()]

        annotated_count = 0

        for contig_id in contig_ids:
            # Get CDS features for this contig (longest isoforms only)
            cur.execute(
                "SELECT Start, \"End\", Strand, Nucleotide_sequence "
                "FROM Contig_annotation_translated "
                "WHERE Contig_id = ? AND Type = 'CDS' "
                "  AND Nucleotide_sequence IS NOT NULL "
                "  AND (Locus_tag IS NULL OR Longest_isoform = true) "
                "ORDER BY Start",
                (contig_id,)
            )
            cds_rows = cur.fetchall()
            if not cds_rows:
                continue

            # Build sorted list of CDS intervals for binary search
            cds_starts = [r[0] for r in cds_rows]
            cds_ends = [r[1] for r in cds_rows]

            # Get all mismatch rows for this contig
            cur.execute(
                "SELECT rowid, First_position, Sequence "
                "FROM Feature_mismatches_with_codons "
                "WHERE Contig_id = ? AND Sequence IS NOT NULL AND Sequence != ''",
                (contig_id,)
            )
            mismatch_rows = cur.fetchall()

            updates = []
            intergenic_updates = []
            for rowid, pos, mismatch_base in mismatch_rows:
                # Find covering CDS via binary search
                # Find rightmost CDS with Start <= pos
                idx = bisect.bisect_right(cds_starts, pos) - 1

                # Check this and any preceding overlapping CDS.
                # Do NOT break early: a CDS with an earlier start can still end after pos
                # even when a later-starting CDS ends before pos (overlapping/nested CDSes).
                cds_found = None
                if idx >= 0:
                    for i in range(idx, -1, -1):
                        if cds_ends[i] >= pos:
                            cds_found = cds_rows[i]
                            break
                        # cds_ends[i] < pos: this CDS ends before pos, but keep searching
                        # backwards — an earlier CDS may be longer and still cover pos

                if cds_found is None:
                    # Position is not inside any CDS — mark as Intergenic
                    intergenic_updates.append(("Intergenic", None, None, rowid))
                    continue

                cds_start, cds_end, strand, nuc_seq = cds_found
                strand = int(strand) if strand is not None else 1

                # Compute offset within CDS
                if strand == -1:
                    offset = cds_end - pos
                    alt_base = complement.get(mismatch_base.upper(), mismatch_base.upper())
                else:
                    offset = pos - cds_start
                    alt_base = mismatch_base.upper()

                codon_idx = offset // 3
                pos_in_codon = offset % 3

                # Extract reference codon from the nucleotide sequence
                codon_start = codon_idx * 3
                codon_end = codon_start + 3
                if codon_end > len(nuc_seq):
                    continue  # Partial codon at end

                ref_codon = nuc_seq[codon_start:codon_end].upper()
                if len(ref_codon) != 3:
                    continue

                # Build mutant codon
                mut_codon = list(ref_codon)
                mut_codon[pos_in_codon] = alt_base
                mut_codon = ''.join(mut_codon)

                # Look up amino acids
                ref_aa = codon_to_aa.get(ref_codon)
                mut_aa = codon_to_aa.get(mut_codon)
                if ref_aa is None or mut_aa is None:
                    continue

                ref_name = codon_to_name.get(ref_codon, ref_aa)
                mut_name = codon_to_name.get(mut_codon, mut_aa)

                category = "Synonymous" if ref_aa == mut_aa else "Non-synonymous"
                codon_change = mut_codon
                aa_change = f"{mut_aa} ({mut_name})"

                updates.append((category, codon_change, aa_change, rowid))
                annotated_count += 1

            # Batch update CDS-overlapping mismatches
            if updates:
                cur.executemany(
                    "UPDATE Feature_mismatches_with_codons "
                    "SET Codon_category = ?, Codon_change = ?, AA_change = ? "
                    "WHERE rowid = ?",
                    updates
                )

            # Batch update intergenic mismatches
            if intergenic_updates:
                cur.executemany(
                    "UPDATE Feature_mismatches_with_codons "
                    "SET Codon_category = ?, Codon_change = ?, AA_change = ? "
                    "WHERE rowid = ?",
                    intergenic_updates
                )

        # Final pass: any mismatch with a non-empty Sequence that still has NULL
        # Codon_category is on a contig/region with no CDS at all → Intergenic
        cur.execute(
            "UPDATE Feature_mismatches_with_codons "
            "SET Codon_category = 'Intergenic' "
            "WHERE Codon_category IS NULL AND Sequence IS NOT NULL AND Sequence != ''"
        )

        db_conn.close()
        print(f"  Annotated {annotated_count} mismatches with codon impact.", flush=True)

    except Exception as e:
        print(f"  WARNING: Could not create mismatch codon annotations: {e}", flush=True)


def calculating_all_features_parallel(list_modules, bam_files, output_db, min_aligned_fraction, min_coverage_depth, curve_ratio, bar_ratio, contig_variation_percentage=0.1, circular=False, n_sample_cores=None, sequencing_type=None, genbank_path=None, assembly_path=None):
    """Process all BAM files in parallel using Rust bindings."""
    if not HAS_RUST:
        sys.exit("ERROR: Rust bindings (thebigbam_rs) are required but not available. Please install them first.")

    if n_sample_cores is None:
        n_sample_cores = max(1, cpu_count() - 1)

    print(f"Using Rust bindings to process {len(bam_files)} samples with rayon ({n_sample_cores} threads)...", flush=True)

    try:
        result = _rust.process_all_samples(
            genbank_path=genbank_path if genbank_path else "",
            bam_files=bam_files,
            output_db=output_db,
            modules=list_modules,
            threads=n_sample_cores,
            sequencing_type=sequencing_type,
            min_aligned_fraction=float(min_aligned_fraction),
            min_coverage_depth=float(min_coverage_depth),
            curve_ratio=float(curve_ratio),
            bar_ratio=float(bar_ratio),
            contig_variation_percentage=float(contig_variation_percentage),
            circular=circular,
            create_indexes=True,
            assembly_path=assembly_path if assembly_path else "",
        )
    except Exception as e:
        print(f"ERROR: Rust processing failed: {e}", flush=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    samples_failed = result.get("samples_failed", 0)
    if samples_failed > 0:
        print(f"Warning: {samples_failed} samples failed to process", flush=True)

def add_calculate_args(parser):
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads available (default: 4)')
    parser.add_argument("-g", "--genbank", help="Path to annotation file: GenBank (.gbk, .gbff) or GFF3 (.gff) format. Required if no BAM files provided.")
    parser.add_argument("-b", "--bam_files", help="Path to bam file or directory containing mapping files (BAM format). Optional if genbank is provided.")
    parser.add_argument("--circular", action="store_true", help="Set if assembly was doubled during mapping (enables modulo logic)")
    parser.add_argument("-o", "--output", required=True, help="Output database file path (.db)")
    parser.add_argument("-m", "--modules", required=False, default=None, help="List of modules to compute (comma-separated). If not provided, all modules are computed. Options: Coverage, Misalignment, Long-reads, Paired-reads, Phage termini")
    parser.add_argument("-a", "--assembly", help="Path to assembly FASTA file (only needed for autoblast when genbank lacks sequence data)")
    parser.add_argument('-s', '--sequencing_type', choices=['long', 'paired-short', 'single-short'], help='Sequencing type (long or short allowed)')
    parser.add_argument("--min_aligned_fraction", type=int, default=50, help="Minimum alignment-length coverage proportion for contig inclusion (default: 50%%)")
    parser.add_argument("--min_coverage_depth", type=int, default=0, help="Minimum mean coverage depth for contig inclusion (disabled by default, e.g. 5 to filter low-depth contigs)")
    parser.add_argument('--variation_percentage', type=float, default=50, help='Run-length encoding ratio for independent features like coverage (default: 50%%)')
    parser.add_argument('--coverage_percentage', type=float, default=10, help='Compressing ratio for features depending on coverage: only values above this %% of the local coverage are kept (default: 10%%)')
    parser.add_argument('--contig_variation_percentage', type=float, default=10, help='Run-length encoding ratio for contig-level features like GC content (default: 10%%)')

VALID_MODULES = ["Coverage", "Misalignment", "Long-reads", "Paired-reads", "Phage termini"]

def run_calculate_args(args):
    genbank_path = getattr(args, 'genbank', None)
    assembly_path = getattr(args, 'assembly', None)

    # Handle optional bam_files
    bam_files = []
    if args.bam_files:
        if os.path.isdir(args.bam_files):
            bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
            if not bam_files:
                print(f"WARNING: No .bam files found in directory '{args.bam_files}'", flush=True)
        elif os.path.isfile(args.bam_files):
            bam_files = [args.bam_files]
        else:
            sys.exit(f"ERROR: BAM path not found: {args.bam_files}")
    
    # Validate: need at least genbank OR bam_files
    if not bam_files and not genbank_path:
        sys.exit("ERROR: You must provide either --bam_files or --genbank (or both).")
    
    if not bam_files:
        print("No BAM files provided. Will only populate contig-level data from GenBank.", flush=True)

    # Handle modules - if no BAM files, modules are ignored
    if args.modules is None:
        # Default: all modules (but they require BAM files)
        requested_modules = VALID_MODULES.copy() if bam_files else []
    else:
        requested_modules = [m.strip() for m in args.modules.split(",")]
        # Validate module names
        for module in requested_modules:
            if module not in VALID_MODULES:
                sys.exit(f"ERROR: Unknown module '{module}'. Valid modules are: {', '.join(VALID_MODULES)}")
        # Warn if modules specified but no BAM files
        if not bam_files and requested_modules:
            print(f"WARNING: Modules {requested_modules} require BAM files - they will be skipped.", flush=True)
            requested_modules = []
    min_aligned_fraction = args.min_aligned_fraction
    min_coverage_depth = args.min_coverage_depth
    variation_percentage = args.variation_percentage
    coverage_percentage = args.coverage_percentage
    contig_variation_percentage = args.contig_variation_percentage
    circular = args.circular
    n_cores = int(args.threads)

    output_db = args.output
    if os.path.exists(output_db):
        sys.exit(f"ERROR: Output file '{output_db}' already exists. Please provide a new path to avoid overwriting.")

    if genbank_path and not os.path.exists(genbank_path):
        sys.exit(f"ERROR: Annotation file not found: {genbank_path}")

    # Validate annotation file extension
    if genbank_path:
        valid_extensions = ('.gbk', '.gbff', '.gb', '.genbank', '.gff', '.gff3')
        if not genbank_path.lower().endswith(valid_extensions):
            sys.exit(f"ERROR: Unsupported annotation file format. Supported extensions: {', '.join(valid_extensions)}")

    if assembly_path and not os.path.exists(assembly_path):
        sys.exit(f"ERROR: Assembly file not found: {assembly_path}")

    print("\nCalculating values for all requested features from mapping files...", flush=True)
    calculating_all_features_parallel(
        requested_modules, bam_files, output_db, min_aligned_fraction, min_coverage_depth, variation_percentage, coverage_percentage,
        contig_variation_percentage=contig_variation_percentage, circular=circular, n_sample_cores=n_cores,
        sequencing_type=args.sequencing_type, genbank_path=genbank_path,
        assembly_path=assembly_path,
    )

    # Store contig sequences for sequence visualization (if source files available)
    if genbank_path or assembly_path:
        print("Storing contig sequences for visualization...", flush=True)
        _store_contig_sequences(output_db, assembly_path=assembly_path, genbank_path=genbank_path)

    # Create translated annotations for CDS visualization
    if genbank_path:
        print("Creating translated annotations for CDS visualization...", flush=True)
        _create_codon_table_and_translated_view(output_db)
        print("Annotating mismatches with codon impact...", flush=True)
        _create_mismatches_with_codons(output_db)

def main():
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    add_calculate_args(parser)
    args = parser.parse_args()
    run_calculate_args(args)


if __name__ == "__main__":
    main()