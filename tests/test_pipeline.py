"""
Integration tests for MGFeatureViewer pipeline.

Tests the calculate command with generated BAM files.
"""

import os
import subprocess
import pytest
import duckdb

# Output database paths
TEST_DB_LINEAR = os.path.join(os.path.dirname(__file__), "test_10kbp_linear.db")
TEST_DB_CIRCULAR = os.path.join(os.path.dirname(__file__), "test_10kbp_circular.db")

def cleanup_dbs():
    """Remove test databases before test run."""
    for db in [TEST_DB_LINEAR, TEST_DB_CIRCULAR]:
        if os.path.exists(db):
            os.remove(db)

def verify_database(db_path: str):
    """Verify database structure and content."""
    assert os.path.exists(db_path), f"Database was not created: {db_path}"

    conn = duckdb.connect(db_path, read_only=True)

    # Check that main tables exist
    tables = conn.execute("SHOW TABLES").fetchall()
    table_names = [t[0] for t in tables]

    assert "Contig" in table_names, "Missing 'Contig' table"
    assert "Sample" in table_names, "Missing 'Sample' table"
    assert "Variable" in table_names, "Missing 'Variable' table"

    # Check that samples were processed
    sample_count = conn.execute("SELECT COUNT(*) FROM Sample").fetchone()[0]
    assert sample_count > 0, "No samples in database"
    print(f"  Processed {sample_count} samples")

    # Check that contigs were processed
    contig_count = conn.execute("SELECT COUNT(*) FROM Contig").fetchone()[0]
    assert contig_count > 0, "No contigs in database"
    print(f"  Found {contig_count} contigs")

    # Check that variables were created
    var_count = conn.execute("SELECT COUNT(*) FROM Variable").fetchone()[0]
    assert var_count > 0, "No variables in database"
    print(f"  Created {var_count} variables")

    # List variables
    variables = conn.execute("SELECT Feature_table_name FROM Variable").fetchall()
    print(f"  Variables: {[v[0] for v in variables]}")

    # Check for expected coverage variable
    var_names = [v[0] for v in variables]
    assert "Feature_primary_reads" in var_names, f"Missing primary_reads variable. Found: {var_names}"

    conn.close()

def run_calculate(bam_dir: str, output_db: str, circular: bool):
    """Run calculate command on BAM files."""
    cmd = [
        "mgfeatureviewer", "calculate",
        "-b", bam_dir,
        "-m", "coverage,phagetermini,assemblycheck",
        "-o", output_db,
        "-t", "4",
        "--coverage_percentage", "1",
        "--variation_percentage", "1",
    ]
    if circular:
        cmd.append("--circular")

    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        print(f"STDOUT:\n{result.stdout}")
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")

    return result

def test_calculate_linear_bams(test_bams, tests_dir):
    """
    Test running calculate on linear (non-circular) BAM files.
    """
    # Clean up before test
    if os.path.exists(TEST_DB_LINEAR):
        os.remove(TEST_DB_LINEAR)

    linear_bams = test_bams["linear"]
    assert len(linear_bams) > 0, "No linear BAM files available"
    print(f"\nTesting with {len(linear_bams)} linear BAMs")

    linear_dir = os.path.join(tests_dir, "linear_bams")
    result = run_calculate(linear_dir, TEST_DB_LINEAR, circular=False)
    assert result.returncode == 0, f"Calculate failed: {result.stderr}"

    print(f"\nVerifying linear database: {TEST_DB_LINEAR}")
    verify_database(TEST_DB_LINEAR)

def test_calculate_circular_bams(test_bams, tests_dir):
    """
    Test running calculate on circular BAM files.
    """
    # Clean up before test
    if os.path.exists(TEST_DB_CIRCULAR):
        os.remove(TEST_DB_CIRCULAR)

    circular_bams = test_bams["circular"]
    assert len(circular_bams) > 0, "No circular BAM files available"
    print(f"\nTesting with {len(circular_bams)} circular BAMs")

    circular_dir = os.path.join(tests_dir, "circular_bams")
    result = run_calculate(circular_dir, TEST_DB_CIRCULAR, circular=True)
    assert result.returncode == 0, f"Calculate failed: {result.stderr}"

    print(f"\nVerifying circular database: {TEST_DB_CIRCULAR}")
    verify_database(TEST_DB_CIRCULAR)
