"""Data accessor module for MGFeatureViewer - pure SQLite."""

import sqlite3
from pathlib import Path


class DataAccessor:
    """Unified data accessor for MGFeatureViewer outputs.

    Supports both:
    - Directory format: directory containing metadata.db
    - Single file format: single .db file

    All feature data is stored in Feature_* tables within the SQLite database.
    """

    def __init__(self, db_path: str):
        """Initialize the data accessor.

        Args:
            db_path: Path to either a directory (containing metadata.db) or .db file
        """
        self.db_path = Path(db_path)

        if self.db_path.is_dir():
            # Directory format: look for metadata.db inside
            self.sqlite_path = self.db_path / "metadata.db"
            if not self.sqlite_path.exists():
                raise FileNotFoundError(f"metadata.db not found in {self.db_path}")
        else:
            # Single file format
            self.sqlite_path = self.db_path
            if not self.sqlite_path.exists():
                raise FileNotFoundError(f"Database file not found: {self.sqlite_path}")

        try:
            self.sqlite_conn = sqlite3.connect(str(self.sqlite_path))
            # Disable WAL mode for compatibility with WSL/Windows cross-filesystem access
            # WAL mode creates -wal and -shm files that don't work well on mounted filesystems
            cur = self.sqlite_conn.cursor()
            cur.execute("PRAGMA journal_mode=DELETE")
            cur.execute("PRAGMA synchronous=NORMAL")
            # Test connection with a simple query
            cur.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1")
            cur.fetchone()
        except sqlite3.Error as e:
            raise RuntimeError(f"Failed to open database {self.sqlite_path}: {e}")

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
            contig_name: Unused, kept for API compatibility
            sample_name: Unused, kept for API compatibility

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

            # Query Feature_* table directly (RLE format: First_position, Last_position, Value)
            cur.execute(
                f"SELECT First_position, Last_position, Value FROM {feature_table} "
                "WHERE Sample_id=? AND Contig_id=? ORDER BY First_position",
                (sample_id, contig_id)
            )
            data_rows = cur.fetchall()
            
            # Debug: log if we have data but it's not being plotted
            if data_rows:
                print(f"[DataAccessor] Found {len(data_rows)} rows for {feature_table} (sample_id={sample_id}, contig_id={contig_id})", flush=True)
                print(f"[DataAccessor] First row: first_pos={data_rows[0][0]}, last_pos={data_rows[0][1]}, value={data_rows[0][2]}", flush=True)
            
            # Expand RLE runs into individual points for plotting
            x_coords = []
            y_coords = []
            for first_pos, last_pos, value in data_rows:
                if type_picked == "bars":
                    # For bars: expand to all positions in the run
                    # (vbar draws one bar per x-coordinate)
                    for pos in range(first_pos, last_pos + 1):
                        x_coords.append(pos)
                        y_coords.append(value)
                else:
                    # For curves: only need start and end points
                    # (line rendering will connect them)
                    if first_pos == last_pos:
                        x_coords.append(first_pos)
                        y_coords.append(value)
                    else:
                        x_coords.extend([first_pos, last_pos])
                        y_coords.extend([value, value])
            
            print(f"[DataAccessor] After expansion: {len(x_coords)} points for {feature_table}", flush=True)
            
            feature_dict["x"] = x_coords
            feature_dict["y"] = y_coords

            # Only append if we have actual data points
            if x_coords:
                list_feature_dict.append(feature_dict)

        return list_feature_dict

    def close(self):
        """Close database connections."""
        if self.sqlite_conn:
            self.sqlite_conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
