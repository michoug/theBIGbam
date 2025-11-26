"""Data accessor module supporting both SQLite (Python output) and DuckDB/Parquet (Rust output)."""

import os
import sqlite3
from pathlib import Path

try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False


class DataAccessor:
    """Unified data accessor for MGFeatureViewer outputs.

    Supports both:
    - Directory format (Rust output): metadata.db + features/*.parquet
    - Single file format (Python output): single .db file with Feature_* tables
    """

    def __init__(self, db_path: str):
        """Initialize the data accessor.

        Args:
            db_path: Path to either a directory (Rust output) or .db file (Python output)
        """
        self.db_path = Path(db_path)

        if self.db_path.is_dir():
            # Rust output: directory with metadata.db and parquet files
            self.mode = "parquet"
            self.sqlite_path = self.db_path / "metadata.db"
            self.features_dir = self.db_path / "features"
            if not self.sqlite_path.exists():
                raise FileNotFoundError(f"metadata.db not found in {self.db_path}")
            if not HAS_DUCKDB:
                raise ImportError("DuckDB is required to read Parquet files. Install with: pip install duckdb")
            self.duckdb_conn = duckdb.connect()
        else:
            # Python output: single SQLite file
            self.mode = "sqlite"
            self.sqlite_path = self.db_path
            self.duckdb_conn = None

        self.sqlite_conn = sqlite3.connect(str(self.sqlite_path))

    def get_sqlite_connection(self) -> sqlite3.Connection:
        """Get the SQLite connection for metadata queries."""
        return self.sqlite_conn

    def get_feature_data(self, feature: str, contig_id: int, sample_id: int,
                         contig_name: str = None, sample_name: str = None):
        """Get feature data for a specific contig and sample.

        Args:
            feature: Feature subplot name (e.g., "Coverage", "Read lengths")
            contig_id: Contig ID from database
            sample_id: Sample ID from database
            contig_name: Contig name (required for parquet mode)
            sample_name: Sample name (required for parquet mode)

        Returns:
            List of feature dicts with x, y data and rendering info
        """
        cur = self.sqlite_conn.cursor()

        # Get rendering info from Variable table
        cur.execute(
            "SELECT Type, Color, Alpha, Fill_alpha, Size, Title, Feature_table_name "
            "FROM Variable WHERE Subplot=?",
            (feature,)
        )
        rows = cur.fetchall()

        list_feature_dict = []
        for row in rows:
            type_picked, color, alpha, fill_alpha, size, title, feature_table = row
            feature_dict = {
                "type": type_picked,
                "color": color,
                "alpha": alpha,
                "fill_alpha": fill_alpha,
                "size": size,
                "title": title,
                "x": [],
                "y": []
            }

            if self.mode == "sqlite":
                # Original SQLite mode - query Feature_* table directly
                cur.execute(
                    f"SELECT Position, Value FROM {feature_table} "
                    "WHERE Sample_id=? AND Contig_id=? ORDER BY Position",
                    (sample_id, contig_id)
                )
                data_rows = cur.fetchall()
                feature_dict["x"] = [r[0] for r in data_rows]
                feature_dict["y"] = [r[1] for r in data_rows]
            else:
                # Parquet mode - query via DuckDB
                if sample_name is None or contig_name is None:
                    raise ValueError("sample_name and contig_name required for parquet mode")

                parquet_file = self.features_dir / f"{sample_name}.parquet"
                if not parquet_file.exists():
                    # No data for this sample
                    list_feature_dict.append(feature_dict)
                    continue

                # Convert Feature_table_name to feature name (strip "Feature_" prefix)
                parquet_feature = feature_table.replace("Feature_", "")

                query = f"""
                    SELECT position, value
                    FROM '{parquet_file}'
                    WHERE contig_name = ? AND feature = ?
                    ORDER BY position
                """
                try:
                    result = self.duckdb_conn.execute(query, [contig_name, parquet_feature]).fetchall()
                    feature_dict["x"] = [r[0] for r in result]
                    feature_dict["y"] = [r[1] for r in result]
                except Exception:
                    # Feature not found in parquet - return empty
                    pass

            list_feature_dict.append(feature_dict)

        return list_feature_dict

    def close(self):
        """Close database connections."""
        if self.sqlite_conn:
            self.sqlite_conn.close()
        if self.duckdb_conn:
            self.duckdb_conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
