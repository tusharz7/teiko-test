#!/usr/bin/env python3
"""Initialise the SQLite database and load cell-count.csv into it.

Usage:
    python load_data.py

Creates (or rebuilds) cell_counts.db next to this script.
"""

import csv
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS cell_populations;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

-- One row per patient. Everything here is constant across a patient's samples.
CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT CHECK (sex IN ('M', 'F')),
    treatment  TEXT NOT NULL,
    response   TEXT CHECK (response IN ('yes', 'no'))  -- NULL when not applicable
);

-- One row per biological sample (a subject sampled at one timepoint).
CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER
);

CREATE TABLE cell_populations (
    population_id INTEGER PRIMARY KEY,
    name          TEXT NOT NULL UNIQUE
);

-- Long format: one row per (sample, population), so adding a new population
-- is an insert rather than a schema change.
CREATE TABLE cell_counts (
    sample_id     TEXT    NOT NULL REFERENCES samples(sample_id),
    population_id INTEGER NOT NULL REFERENCES cell_populations(population_id),
    count         INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
) WITHOUT ROWID;

CREATE INDEX idx_subjects_project   ON subjects(project_id);
CREATE INDEX idx_subjects_filters   ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_subject    ON samples(subject_id);
CREATE INDEX idx_samples_type_time  ON samples(sample_type, time_from_treatment_start);
CREATE INDEX idx_counts_population  ON cell_counts(population_id);
"""


def to_int(value):
    value = value.strip()
    return int(value) if value else None


def to_text(value):
    value = value.strip()
    return value or None


def load(csv_path=CSV_PATH, db_path=DB_PATH):
    if not csv_path.exists():
        sys.exit(f"error: {csv_path} not found")

    projects = set()
    subjects = {}
    samples = []
    counts = []
    population_ids = {name: i for i, name in enumerate(POPULATIONS, start=1)}

    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            projects.add(row["project"])

            subject = (
                row["subject"],
                row["project"],
                row["condition"].strip(),
                to_int(row["age"]),
                to_text(row["sex"]),
                row["treatment"].strip(),
                to_text(row["response"]),
            )
            previous = subjects.setdefault(row["subject"], subject)
            if previous != subject:
                sys.exit(
                    f"error: subject {row['subject']} has conflicting metadata "
                    f"at sample {row['sample']}"
                )

            samples.append(
                (
                    row["sample"],
                    row["subject"],
                    row["sample_type"].strip(),
                    to_int(row["time_from_treatment_start"]),
                )
            )
            for name in POPULATIONS:
                counts.append((row["sample"], population_ids[name], int(row[name])))

    conn = sqlite3.connect(db_path)
    try:
        conn.executescript(SCHEMA)
        with conn:
            conn.executemany(
                "INSERT INTO projects VALUES (?)", [(p,) for p in sorted(projects)]
            )
            conn.executemany(
                "INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)", subjects.values()
            )
            conn.executemany("INSERT INTO samples VALUES (?, ?, ?, ?)", samples)
            conn.executemany(
                "INSERT INTO cell_populations VALUES (?, ?)",
                [(i, name) for name, i in population_ids.items()],
            )
            conn.executemany("INSERT INTO cell_counts VALUES (?, ?, ?)", counts)
    finally:
        conn.close()

    print(
        f"Loaded {len(projects)} projects, {len(subjects)} subjects, "
        f"{len(samples)} samples, {len(counts)} cell counts into {db_path.name}"
    )


if __name__ == "__main__":
    load()
