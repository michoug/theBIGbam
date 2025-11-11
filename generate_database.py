import sqlite3
import sys

if len(sys.argv) != 2:
    print("Usage: python generate_database.py <database_name.db>")
    sys.exit(1)

DB = sys.argv[1]
conn = sqlite3.connect(DB)
cur = conn.cursor()

# --- Drop old tables if they exist ---
cur.executescript("""
DROP TABLE IF EXISTS Contig;
DROP TABLE IF EXISTS Sample;
DROP TABLE IF EXISTS Variable;
""")

# --- Create Contig and Sample tables ---
cur.execute("""
CREATE TABLE Contig (
    Contig_id INTEGER PRIMARY KEY AUTOINCREMENT,
    Contig_name TEXT UNIQUE,
    Contig_length INTEGER
);
""")

cur.execute("""
CREATE TABLE Sample (
    Sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
    Sample_name TEXT UNIQUE
);
""")

# --- Create Variable table ---
cur.execute("""
CREATE TABLE Variable (
    Variable_id INTEGER PRIMARY KEY AUTOINCREMENT,
    Variable_name TEXT,
    Status TEXT,
    Type TEXT,
    Color TEXT,
    Size REAL,
    Feature_table_name TEXT
);
""")

conn.commit()
conn.close()
print(f"Empty database created: {DB}")