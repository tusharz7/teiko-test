#!/usr/bin/env python3
"""Single-file analysis module for Bob Loblaw's immune cell population study.

Covers all four parts of the assignment end to end:
    Part 1 - database schema + CSV loading
    Part 2 - per-sample population frequencies
    Part 3 - responder vs. non-responder comparison (melanoma / miraclib / PBMC)
    Part 4 - baseline (day 0) subset breakdown

Every query lives behind a small function that returns a pandas DataFrame, so
this module doubles as the data layer for dashboard.py and as a standalone
CLI report (`python analysis.py`). The authoritative loader required by the
assignment spec is still load_data.py; build_database() below is the same
logic, kept here so the whole pipeline can be read/run from one file.
"""

import csv
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"

# The five measured immune cell populations, in a fixed display order used
# throughout (table ordering, boxplot x-axis, stats output).
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


# The cohort the assignment asks about in Parts 3 and 4. The query functions
# take these as defaults so the dashboard can point them at other cohorts.
DEFAULT_CONDITION = "melanoma"
DEFAULT_TREATMENT = "miraclib"
DEFAULT_SAMPLE_TYPE = "PBMC"


def connect(db_path=DB_PATH):
    """Open a connection to the cell-count database."""
    return sqlite3.connect(db_path)


# =============================================================================
# Part 1: Data Management - schema definition and CSV loading
# =============================================================================
#
# Normalized into four entities so each fact is stored once:
#   projects   - top-level trial grouping
#   subjects   - one row per patient; demographics/treatment are constant
#                across all of a patient's samples, so they don't belong on
#                the sample or count rows
#   samples    - one row per (subject, timepoint) biological sample
#   cell_counts - long format (sample, population) -> count, rather than one
#                 column per population; adding a new population later is an
#                 insert, not a schema migration

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


def _to_int(value):
    """Parse an optional integer field; blank CSV cells become NULL."""
    value = value.strip()
    return int(value) if value else None


def _to_text(value):
    """Parse an optional text field; blank CSV cells become NULL."""
    value = value.strip()
    return value or None


def build_database(csv_path=CSV_PATH, db_path=DB_PATH):
    """Create (or rebuild) the SQLite database and load cell-count.csv into it.

    Subject-level fields (condition, age, sex, treatment, response) are
    de-duplicated across a subject's rows and checked for consistency, since
    the source CSV repeats them once per sample.
    """
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
                _to_int(row["age"]),
                _to_text(row["sex"]),
                row["treatment"].strip(),
                _to_text(row["response"]),
            )
            # Guard against a malformed CSV where the same subject appears
            # with conflicting demographics across different sample rows.
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
                    _to_int(row["time_from_treatment_start"]),
                )
            )
            for name in POPULATIONS:
                counts.append((row["sample"], population_ids[name], int(row[name])))

    conn = sqlite3.connect(db_path)
    try:
        conn.executescript(SCHEMA)
        with conn:  # commits on success, rolls back on exception
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

    return {
        "projects": len(projects),
        "subjects": len(subjects),
        "samples": len(samples),
        "counts": len(counts),
    }


# =============================================================================
# Part 2: Initial Analysis - relative frequency of each population per sample
# =============================================================================

FREQUENCIES_QUERY = """
SELECT
    cc.sample_id AS sample,
    totals.total_count AS total_count,
    cp.name AS population,
    cc.count AS count,
    ROUND(100.0 * cc.count / totals.total_count, 2) AS percentage
FROM cell_counts cc
JOIN cell_populations cp ON cp.population_id = cc.population_id
JOIN (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
) totals ON totals.sample_id = cc.sample_id
ORDER BY cc.sample_id, cp.population_id;
"""


def frequencies(db_path=DB_PATH):
    """One row per (sample, population): count and % of that sample's total."""
    with connect(db_path) as conn:
        return pd.read_sql_query(FREQUENCIES_QUERY, conn)


# =============================================================================
# Part 3: Statistical Analysis - responders vs. non-responders
#         (melanoma, miraclib, PBMC samples only)
# =============================================================================

RESPONDER_QUERY = """
SELECT
    cc.sample_id AS sample,
    sa.subject_id AS subject,
    su.response AS response,
    totals.total_count AS total_count,
    cp.name AS population,
    cc.count AS count,
    ROUND(100.0 * cc.count / totals.total_count, 2) AS percentage
FROM cell_counts cc
JOIN cell_populations cp ON cp.population_id = cc.population_id
JOIN samples sa ON sa.sample_id = cc.sample_id
JOIN subjects su ON su.subject_id = sa.subject_id
JOIN (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
) totals ON totals.sample_id = cc.sample_id
WHERE su.condition = ?
  AND su.treatment = ?
  AND sa.sample_type = ?
  AND su.response IN ('yes', 'no')
ORDER BY cp.population_id, cc.sample_id;
"""


def responder_frequencies(
    db_path=DB_PATH,
    condition=DEFAULT_CONDITION,
    treatment=DEFAULT_TREATMENT,
    sample_type=DEFAULT_SAMPLE_TYPE,
):
    """Long-format population frequencies for one cohort's samples (by default
    melanoma/miraclib/PBMC), labelled by responder status. Feeds both the
    boxplot and the stats test.
    """
    with connect(db_path) as conn:
        return pd.read_sql_query(
            RESPONDER_QUERY, conn, params=(condition, treatment, sample_type)
        )


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR-adjusted p-values (q-values).

    Implemented directly (rather than pulling in statsmodels) since it's a
    handful of lines: rank p-values ascending, scale each by n/rank, then
    enforce monotonicity by taking a running minimum from the largest
    p-value down.
    """
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out


def responder_stats(df=None, db_path=DB_PATH):
    """Per-population Mann-Whitney U test, responders vs. non-responders.

    Mann-Whitney (rather than a t-test) because relative frequencies are
    bounded percentages with no guarantee of normality, and it's robust to
    the unequal group sizes typical of a responder/non-responder split.
    Testing all 5 populations simultaneously inflates false positives, so
    p-values are Benjamini-Hochberg corrected before flagging significance.
    """
    if df is None:
        df = responder_frequencies(db_path)

    rows = []
    for population in POPULATIONS:
        subset = df[df["population"] == population]
        responders = subset.loc[subset["response"] == "yes", "percentage"]
        non_responders = subset.loc[subset["response"] == "no", "percentage"]
        if responders.empty or non_responders.empty:
            # Nothing to compare (e.g. a cohort with only one response group).
            continue
        stat, p_value = mannwhitneyu(
            responders, non_responders, alternative="two-sided"
        )
        rows.append(
            {
                "population": population,
                "n_responders": len(responders),
                "n_non_responders": len(non_responders),
                "median_responders": responders.median(),
                "median_non_responders": non_responders.median(),
                "u_stat": stat,
                "p_value": p_value,
            }
        )

    result = pd.DataFrame(
        rows,
        columns=[
            "population",
            "n_responders",
            "n_non_responders",
            "median_responders",
            "median_non_responders",
            "u_stat",
            "p_value",
        ],
    )
    result["p_adj"] = benjamini_hochberg(result["p_value"].to_numpy())
    result["significant"] = result["p_adj"] < 0.05
    return result.sort_values("p_adj").reset_index(drop=True)


# =============================================================================
# Part 4: Data Subset Analysis - baseline (day 0) melanoma/PBMC/miraclib subset
# =============================================================================

BASELINE_QUERY = """
SELECT
    sa.sample_id AS sample,
    sa.subject_id AS subject,
    su.project_id AS project,
    su.response AS response,
    su.sex AS sex
FROM samples sa
JOIN subjects su ON su.subject_id = sa.subject_id
WHERE su.condition = ?
  AND su.treatment = ?
  AND sa.sample_type = ?
  AND sa.time_from_treatment_start = 0
ORDER BY sa.sample_id;
"""


def baseline_subset(
    db_path=DB_PATH,
    condition=DEFAULT_CONDITION,
    treatment=DEFAULT_TREATMENT,
    sample_type=DEFAULT_SAMPLE_TYPE,
):
    """One cohort's samples taken at baseline (time 0); by default
    melanoma / miraclib / PBMC.
    """
    with connect(db_path) as conn:
        return pd.read_sql_query(
            BASELINE_QUERY, conn, params=(condition, treatment, sample_type)
        )


def baseline_summary(df=None, db_path=DB_PATH):
    """Breakdown of the baseline subset: samples per project, and subject
    counts by response and by sex.

    Response/sex are subject-level attributes, so subjects are de-duplicated
    first - otherwise a subject with two qualifying baseline samples would
    be double-counted.
    """
    if df is None:
        df = baseline_subset(db_path)

    samples_per_project = (
        df.groupby("project")["sample"].nunique().reset_index(name="n_samples")
    )

    subjects = df.drop_duplicates(subset="subject")
    response_counts = (
        subjects.groupby("response")["subject"]
        .nunique()
        .reset_index(name="n_subjects")
    )
    sex_counts = (
        subjects.groupby("sex")["subject"].nunique().reset_index(name="n_subjects")
    )

    return samples_per_project, response_counts, sex_counts


# =============================================================================
# Dashboard helpers - dataset-wide overview and the available cohort filters
# =============================================================================

COHORT_QUERY = """
SELECT
    su.project_id AS project,
    su.condition AS condition,
    su.treatment AS treatment,
    sa.sample_type AS sample_type,
    COUNT(DISTINCT su.subject_id) AS n_subjects,
    COUNT(*) AS n_samples
FROM samples sa
JOIN subjects su ON su.subject_id = sa.subject_id
GROUP BY su.project_id, su.condition, su.treatment, sa.sample_type
ORDER BY su.project_id, su.condition, su.treatment, sa.sample_type;
"""


def cohort_summary(db_path=DB_PATH):
    """Subject and sample counts per project / condition / treatment /
    sample type. Drives the overview tab and the cohort filter options.
    """
    with connect(db_path) as conn:
        return pd.read_sql_query(COHORT_QUERY, conn)


def dataset_totals(db_path=DB_PATH):
    """Headline counts for the whole database."""
    with connect(db_path) as conn:
        row = conn.execute(
            """
            SELECT
                (SELECT COUNT(*) FROM projects),
                (SELECT COUNT(*) FROM subjects),
                (SELECT COUNT(*) FROM samples),
                (SELECT SUM(count) FROM cell_counts)
            """
        ).fetchone()
    return dict(zip(["projects", "subjects", "samples", "cells"], row))


# =============================================================================
# CLI entry point - runs the full pipeline and prints a report for each part
# =============================================================================

if __name__ == "__main__":
    # load_data.py is the loader; only build here if it hasn't been run.
    if not DB_PATH.exists():
        stats = build_database()
        print(
            f"Loaded {stats['projects']} projects, {stats['subjects']} subjects, "
            f"{stats['samples']} samples, {stats['counts']} cell counts "
            f"into {DB_PATH.name}"
        )

    OUTPUT_DIR.mkdir(exist_ok=True)

    def save(df, name):
        df.to_csv(OUTPUT_DIR / name, index=False)
        print(f"  wrote {OUTPUT_DIR.name}/{name} ({len(df)} rows)")

    print("\nPart 2: frequencies")
    freq = frequencies()
    print(freq.head())
    save(freq, "part2_frequencies.csv")

    print("\nPart 3: responder stats")
    responders = responder_frequencies()
    stats_table = responder_stats(responders)
    print(stats_table)
    save(responders, "part3_responder_frequencies.csv")
    save(stats_table, "part3_responder_stats.csv")

    print("\nPart 4: baseline summary")
    baseline = baseline_subset()
    save(baseline, "part4_baseline_samples.csv")
    for name, result in zip(
        ["samples_per_project", "response_counts", "sex_counts"],
        baseline_summary(baseline),
    ):
        print(f"\n{name.replace('_', ' ')}:")
        print(result)
        save(result, f"part4_{name}.csv")
