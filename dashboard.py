#!/usr/bin/env python3
"""Interactive dashboard for Bob's immune cell population analysis.

Usage:
    streamlit run dashboard.py
"""

from pathlib import Path

import plotly.express as px
import streamlit as st

import analysis

ROOT = Path(__file__).resolve().parent

# Palette taken from teiko.bio (page colours are set in .streamlit/config.toml).
RED = "#FF3534"
DEEP_RED = "#C8102E"
INDIGO = "#6C7FD8"
INK = "#1A1A1F"
GREY = "#9494A0"
BORDER = "#E9E7E3"

RESPONSE_LABELS = {"yes": "Responder", "no": "Non-Responder"}
RESPONSE_COLOR_MAP = {"Responder": INDIGO, "Non-Responder": RED}
TREATMENT_COLOR_MAP = {"miraclib": RED, "phauximab": INDIGO, "none": GREY}
SEX_LABELS = {"F": "Female", "M": "Male"}

st.set_page_config(page_title="Loblaw Bio - Cell Population Analysis", layout="wide")

st.markdown(
    f"""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Montserrat:wght@600;700&display=swap');
    html, body, [class*="st-"], .stMarkdown {{ font-family: 'Inter', sans-serif; }}
    /* Icons are drawn with an icon font; keep it, or their names show as text. */
    [data-testid="stIconMaterial"] {{ font-family: 'Material Symbols Rounded' !important; }}
    h1, h2, h3, h4 {{ font-family: 'Montserrat', sans-serif !important; color: {INK}; }}
    [data-testid="stMetric"] {{
        background: #FFFFFF; border: 1px solid {BORDER}; border-radius: 16px;
        padding: 1rem 1.25rem;
    }}
    [data-testid="stMetricLabel"] p {{
        color: {INK}; font-size: 0.75rem; font-weight: 600;
        letter-spacing: 0.06em; text-transform: uppercase;
    }}
    [data-testid="stMetricValue"], [data-testid="stMetricValue"] * {{
        font-family: 'Montserrat', sans-serif; font-weight: 700; color: {INK} !important;
    }}
    [data-testid="stVerticalBlockBorderWrapper"] {{ background: #FFFFFF; border-radius: 16px; }}
    .stTabs [aria-selected="true"] {{ font-weight: 600; }}
    </style>
    """,
    unsafe_allow_html=True,
)


def style(fig, height=380):
    """Apply the shared chart look."""
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=10, t=30, b=10),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#FFFFFF",
        font=dict(family="Inter, sans-serif", color=INK),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, title_text=""),
    )
    fig.update_xaxes(showgrid=False, linecolor=BORDER)
    fig.update_yaxes(gridcolor="#F0EEEA", linecolor=BORDER, zeroline=False)
    return fig


def chart(fig):
    st.plotly_chart(fig, use_container_width=True)


def table(df, **kwargs):
    st.dataframe(df, width="stretch", hide_index=True, **kwargs)


def count_bar(df, x, y, color=INDIGO):
    """Small single-series bar chart with the value printed on each bar."""
    fig = px.bar(df, x=x, y=y, text=y, color_discrete_sequence=[color])
    fig.update_traces(textposition="outside", cliponaxis=False)
    fig.update_xaxes(type="category", title_text="")
    fig.update_yaxes(title_text="")
    return style(fig, height=260)


db_path = ROOT / "cell_counts.db"
if not db_path.exists():
    # `make pipeline` normally creates the database. On a fresh checkout (e.g.
    # the hosted dashboard) build it from the CSV with the same loader logic.
    with st.spinner("Building the database from cell-count.csv..."):
        analysis.build_database(db_path=db_path)

load_frequencies = st.cache_data(analysis.frequencies)
load_cohorts = st.cache_data(analysis.cohort_summary)
load_totals = st.cache_data(analysis.dataset_totals)
load_responders = st.cache_data(analysis.responder_frequencies)
load_baseline = st.cache_data(analysis.baseline_subset)

cohorts = load_cohorts(db_path)
totals = load_totals(db_path)

# Parts 3 and 4 look at the cohort Bob asked about.
condition = analysis.DEFAULT_CONDITION
treatment = analysis.DEFAULT_TREATMENT
sample_type = analysis.DEFAULT_SAMPLE_TYPE
cohort_label = f"{condition} · {treatment} · {sample_type}"

st.title("Immune Cell Population Analysis")

overview_tab, tab2, tab3, tab4 = st.tabs(
    [
        "Overview",
        "Part 2 - Population Frequencies",
        "Part 3 - Responders vs Non-Responders",
        "Part 4 - Baseline Subset",
    ]
)

# ---------------------------------------------------------------------------
# Overview
# ---------------------------------------------------------------------------
with overview_tab:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Projects", f"{totals['projects']:,}")
    c2.metric("Subjects", f"{totals['subjects']:,}")
    c3.metric("Samples", f"{totals['samples']:,}")
    c4.metric("Cells Counted", f"{totals['cells'] / 1e6:,.0f}M")

    left, right = st.columns(2)
    with left, st.container(border=True):
        st.subheader("Samples by Condition and Treatment")
        by_arm = cohorts.groupby(["condition", "treatment"], as_index=False)[
            ["n_subjects", "n_samples"]
        ].sum()
        # Untreated healthy controls go last, after the disease cohorts.
        by_arm = by_arm.sort_values(
            "condition", key=lambda c: c == "healthy", kind="stable"
        )
        table(
            by_arm,
            column_config={
                "condition": "Condition",
                "treatment": "Treatment",
                "n_subjects": "Subjects",
                "n_samples": "Samples",
            },
        )
        fig = px.bar(
            by_arm,
            x="condition",
            y="n_samples",
            color="treatment",
            color_discrete_map=TREATMENT_COLOR_MAP,
            labels={"condition": "", "n_samples": "Samples"},
        )
        chart(style(fig, height=300))

    with right, st.container(border=True):
        st.subheader("Samples by Project and Sample Type")
        by_project = cohorts.groupby(["project", "sample_type"], as_index=False)[
            ["n_subjects", "n_samples"]
        ].sum()
        table(
            by_project,
            column_config={
                "project": "Project",
                "sample_type": "Sample Type",
                "n_subjects": "Subjects",
                "n_samples": "Samples",
            },
        )
        fig = px.bar(
            by_project,
            x="project",
            y="n_samples",
            color="sample_type",
            color_discrete_map={"PBMC": INDIGO, "WB": INK},
            labels={"project": "", "n_samples": "Samples"},
        )
        chart(style(fig, height=300))

# ---------------------------------------------------------------------------
# Part 2
# ---------------------------------------------------------------------------
with tab2:
    st.subheader("Relative Frequency of Each Cell Population per Sample")
    freq_df = load_frequencies(db_path)

    samples = sorted(freq_df["sample"].unique())
    # The filter widget is drawn at the bottom of the tab, so read its value
    # from session state here.
    selected = st.session_state.get("sample_filter", [])
    view = freq_df[freq_df["sample"].isin(selected)] if selected else freq_df

    per_sample = view.drop_duplicates("sample")
    c1, c2, c3 = st.columns(3)
    c1.metric("Samples Shown", f"{len(per_sample):,}")
    c2.metric("Populations", f"{view['population'].nunique()}")
    c3.metric("Median Cells per Sample", f"{per_sample['total_count'].median():,.0f}")

    left, right = st.columns([3, 2])
    with left, st.container(border=True):
        st.markdown("**Frequency Table**")
        table(
            view,
            height=420,
            column_config={
                "sample": "Sample",
                "total_count": "Total Count",
                "population": "Population",
                "count": "Count",
                "percentage": "Percentage",
            },
        )
    with right, st.container(border=True):
        if selected and len(selected) <= 25:
            st.markdown("**Composition of Selected Samples**")
            fig = px.bar(
                view,
                x="percentage",
                y="sample",
                color="population",
                orientation="h",
                color_discrete_sequence=[INDIGO, RED, INK, GREY, "#E8B4A0"],
                category_orders={"population": analysis.POPULATIONS},
                labels={
                    "percentage": "Relative Frequency (%)",
                    "sample": "",
                    "population": "Population",
                },
            )
            fig.update_layout(barmode="stack")
        else:
            st.markdown("**Mean Relative Frequency Across Samples Shown**")
            mean_df = (
                view.groupby("population", as_index=False)["percentage"].mean().round(2)
            )
            fig = px.bar(
                mean_df,
                x="population",
                y="percentage",
                text="percentage",
                color_discrete_sequence=[INDIGO],
                category_orders={"population": analysis.POPULATIONS},
                labels={"population": "", "percentage": "Mean Relative Frequency (%)"},
            )
            fig.update_traces(textposition="outside", cliponaxis=False)
        chart(style(fig, height=420))

    st.multiselect(
        "Filter by sample (leave empty to show all)", samples, key="sample_filter"
    )

# ---------------------------------------------------------------------------
# Part 3
# ---------------------------------------------------------------------------
with tab3:
    st.subheader("Responders vs Non-Responders")
    st.caption(f"Cohort: {cohort_label}")
    resp_df = load_responders(db_path, condition, treatment, sample_type)
    stats_df = analysis.responder_stats(resp_df)

    if stats_df.empty:
        st.warning(
            f"No responder and non-responder samples to compare for {cohort_label}."
        )
    else:
        resp_df = resp_df.assign(response=resp_df["response"].map(RESPONSE_LABELS))
        significant = stats_df[stats_df["significant"]]
        per_sample = resp_df.drop_duplicates("sample")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Samples", f"{len(per_sample):,}")
        c2.metric("Responder Samples", f"{(per_sample['response'] == 'Responder').sum():,}")
        c3.metric(
            "Non-Responder Samples",
            f"{(per_sample['response'] == 'Non-Responder').sum():,}",
        )
        c4.metric("Significant Populations", f"{len(significant)} of {len(stats_df)}")

        if significant.empty:
            st.info(
                "No cell population shows a statistically significant difference "
                "between responders and non-responders after FDR correction "
                "(adjusted p < 0.05)."
            )
        else:
            names = ", ".join(significant["population"])
            st.success(f"Significant difference (adjusted p < 0.05) in: {names}")

        with st.container(border=True):
            st.markdown("**Populations Ranked by Adjusted P-Value**")
            st.caption("Mann-Whitney U test per population, Benjamini-Hochberg adjusted.")
            table(
                stats_df,
                column_config={
                    "population": "Population",
                    "n_responders": "Responder Samples",
                    "n_non_responders": "Non-Responder Samples",
                    "median_responders": st.column_config.NumberColumn(
                        "Median % (Responders)", format="%.2f"
                    ),
                    "median_non_responders": st.column_config.NumberColumn(
                        "Median % (Non-Responders)", format="%.2f"
                    ),
                    "u_stat": st.column_config.NumberColumn("U Statistic", format="%.0f"),
                    "p_value": st.column_config.NumberColumn("P-Value", format="%.4f"),
                    "p_adj": st.column_config.NumberColumn("Adjusted P", format="%.4f"),
                    "significant": "Significant",
                },
            )

        with st.container(border=True):
            st.markdown("**Relative Frequency by Population**")
            fig = px.box(
                resp_df,
                x="population",
                y="percentage",
                color="response",
                color_discrete_map=RESPONSE_COLOR_MAP,
                points="all",
                labels={
                    "population": "",
                    "percentage": "Relative Frequency (%)",
                    "response": "Response",
                },
                category_orders={
                    "population": analysis.POPULATIONS,
                    "response": list(RESPONSE_COLOR_MAP),
                },
            )
            fig.update_traces(marker=dict(size=3, opacity=0.35))
            fig.update_layout(boxmode="group")
            chart(style(fig, height=480))

# ---------------------------------------------------------------------------
# Part 4
# ---------------------------------------------------------------------------
with tab4:
    st.subheader("Baseline (Day 0) Samples")
    st.caption(f"Cohort: {cohort_label}")
    baseline_df = load_baseline(db_path, condition, treatment, sample_type)

    if baseline_df.empty:
        st.warning(
            f"No baseline samples for {cohort_label}."
        )
    else:
        samples_per_project, response_counts, sex_counts = analysis.baseline_summary(
            baseline_df
        )
        response_counts = response_counts.assign(
            response=response_counts["response"].map(RESPONSE_LABELS)
        )
        sex_counts = sex_counts.assign(sex=sex_counts["sex"].map(SEX_LABELS))

        c1, c2, c3 = st.columns(3)
        c1.metric("Baseline Samples", f"{baseline_df['sample'].nunique():,}")
        c2.metric("Subjects", f"{baseline_df['subject'].nunique():,}")
        c3.metric("Projects", f"{baseline_df['project'].nunique():,}")

        col1, col2, col3 = st.columns(3)
        with col1, st.container(border=True):
            st.markdown("**Samples per Project**")
            chart(count_bar(samples_per_project, "project", "n_samples", INK))
            table(
                samples_per_project,
                column_config={"project": "Project", "n_samples": "Samples"},
            )
        with col2, st.container(border=True):
            st.markdown("**Responders vs Non-Responders**")
            if response_counts.empty:
                st.caption("No response recorded for this cohort.")
            else:
                chart(count_bar(response_counts, "response", "n_subjects", INDIGO))
                table(
                    response_counts,
                    column_config={"response": "Response", "n_subjects": "Subjects"},
                )
        with col3, st.container(border=True):
            st.markdown("**Subjects by Sex**")
            chart(count_bar(sex_counts, "sex", "n_subjects", RED))
            table(sex_counts, column_config={"sex": "Sex", "n_subjects": "Subjects"})

        with st.expander("Show the Baseline Samples"):
            table(
                baseline_df,
                column_config={
                    "sample": "Sample",
                    "subject": "Subject",
                    "project": "Project",
                    "response": "Response",
                    "sex": "Sex",
                },
            )
