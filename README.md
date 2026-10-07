# Loblaw Bio: Immune Cell Population Analysis

Analysis of immune cell counts from Bob Loblaw's clinical trial: a SQLite
database built from `cell-count.csv`, a Python analysis of how the drug
candidate affects five immune cell populations, and an interactive dashboard
that presents the results.

**Dashboard:** <!-- TODO: replace with the deployed URL --> _link to be added after deployment_

## Running it

The project runs in GitHub Codespaces (or any machine with Python 3.9+ and
`make`). From the repository root:

```bash
make setup      # install dependencies from requirements.txt
make pipeline   # build cell_counts.db, run the analysis, write outputs/
make dashboard  # start the dashboard on http://localhost:8501
```

In Codespaces, `make dashboard` forwards port 8501; open it from the
notification or the **Ports** tab.

`make pipeline` runs two steps in order, with no manual input:

1. `python load_data.py` creates `cell_counts.db` in the repository root and
   loads every row of `cell-count.csv`.
2. `python analysis.py` runs Parts 2-4 against that database, prints a report,
   and writes the result tables to `outputs/`:

| File | Contents |
|---|---|
| `part2_frequencies.csv` | Relative frequency of each population in each sample |
| `part3_responder_frequencies.csv` | Frequencies for melanoma / miraclib / PBMC samples, labelled by response |
| `part3_responder_stats.csv` | Responder vs non-responder test per population |
| `part4_baseline_samples.csv` | Baseline melanoma / miraclib / PBMC samples |
| `part4_samples_per_project.csv` | Baseline samples per project |
| `part4_response_counts.csv` | Baseline subjects by response |
| `part4_sex_counts.csv` | Baseline subjects by sex |

`cell_counts.db` and `outputs/` are generated, so they are not committed.

## Database schema

```
projects (project_id)
    |
subjects (subject_id, project_id, condition, age, sex, treatment, response)
    |
samples (sample_id, subject_id, sample_type, time_from_treatment_start)
    |
cell_counts (sample_id, population_id, count)  ---  cell_populations (population_id, name)
```

- **`projects`**: one row per project.
- **`subjects`**: one row per patient. Condition, age, sex, treatment and
  response are the same on every one of a patient's samples, so they are stored
  once here. The loader stops with an error if a file ever contradicts that.
  Response is `NULL` for subjects with no recorded response (the untreated
  healthy controls).
- **`samples`**: one row per biological sample, i.e. one patient at one
  timepoint.
- **`cell_populations`** and **`cell_counts`**: counts are stored in long
  format, one row per sample and population.

### Why this design, and how it scales

- **Each fact is stored once.** Correcting a patient's response is a one-row
  update, not an update to every sample row.
- **New populations need no schema change.** Adding a sixth cell population, or
  a few hundred, is an insert into `cell_populations` plus rows in
  `cell_counts`. Analysis queries group by population, so they keep working
  unchanged.
- **Hundreds of projects and thousands of samples.** Projects, subjects and
  samples are separate tables joined on indexed keys, and the columns the
  analyses filter on (condition, treatment, response, sample type, timepoint)
  are indexed, so cohort queries stay index lookups as the data grows.
  `cell_counts` is keyed on `(sample_id, population_id)` and stored without a
  separate row id, which keeps the largest table compact.
- **New kinds of analysis.** Because metadata lives at the level it belongs to,
  questions such as "frequency over time per subject" or "compare treatments
  within a condition" are joins and `GROUP BY`s over the existing tables.
- **Beyond SQLite.** The schema is plain SQL with foreign keys and would move
  to PostgreSQL as is if many people needed to write to it at once. With more
  metadata, treatments and sample types would become their own lookup tables.

## Code structure

| File | Role |
|---|---|
| `load_data.py` | Part 1. Creates the schema and loads the CSV. Standard library only. |
| `analysis.py` | Parts 2-4. Each analysis is a SQL query behind a small function that returns a pandas DataFrame. Run as a script, it prints a report and writes `outputs/`. |
| `dashboard.py` | Streamlit dashboard. Imports the functions in `analysis.py`, so the dashboard and the pipeline show the same numbers. |
| `Makefile` | `setup`, `pipeline` and `dashboard` targets. |
| `.streamlit/config.toml` | Dashboard colour theme. |

The analysis is kept separate from the dashboard so that every number has one
source: the dashboard contains layout and charts, not calculations. Heavy
lifting (totals, percentages, cohort filters) is done in SQL, where it is
closest to the data; pandas and SciPy are used only for the statistics.

## Analysis notes

- **Part 2**: for each sample, the total cell count is the sum of the five
  populations, and each population's relative frequency is its count as a
  percentage of that total.
- **Part 3**: melanoma patients on miraclib, PBMC samples only, responders
  against non-responders. Each population is compared with a two-sided
  Mann-Whitney U test, because relative frequencies are bounded percentages
  that need not be normally distributed. The five p-values are
  Benjamini-Hochberg adjusted, and a population is called significant at
  adjusted p < 0.05.
  **Result:** no population differs significantly after adjustment. CD4 T
  cells come closest (p = 0.013, adjusted p = 0.067), with a slightly higher
  median frequency in responders (30.22% against 29.66%).
- **Part 4**: baseline (time 0) melanoma PBMC samples from miraclib-treated
  patients: 656 samples from 656 subjects, 384 in prj1 and 272 in prj3;
  331 responders and 325 non-responders; 344 male and 312 female.

## Dashboard

Four tabs: an overview of the dataset, then one tab per part.

- **Overview**: totals, and samples by condition, treatment, project and
  sample type.
- **Part 2**: the frequency table, with a sample filter and a composition chart.
- **Part 3**: the statistics table ranked by adjusted p-value, and boxplots of
  each population for responders and non-responders.
- **Part 4**: the baseline subset broken down by project, response and sex.

The dashboard reads `cell_counts.db` as created by `make pipeline`. If the
file is missing, it builds the database from `cell-count.csv` on first load.
