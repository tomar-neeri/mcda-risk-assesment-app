import streamlit as st
import pandas as pd
import numpy as np
import re
import seaborn as sns
import matplotlib.pyplot as plt
import os
from io import BytesIO

st.set_page_config(layout="wide")

# --- Layout: Logo + Title ---
col1, col2 = st.columns([1, 5])
with col1:
    if os.path.exists("logo.png"):
        st.image("logo.png", width=120)
with col2:
    st.markdown("## Graphical User Interface Tool for Assessing and Prioritizing Antimicrobial Resistance Risks from Metagenomic Datasets")

# --- Load README.md instructions ---
def load_readme_md():
    pr = os.path.dirname(__file__)
    readme_path = os.path.join(pr, "README.md")
    if os.path.isfile(readme_path):
        with open(readme_path, "r", encoding="utf-8") as f:
            return f.read()
    return "README.md file not found."

readme_content = load_readme_md()

with st.expander("User Instructions", expanded=False):
    st.markdown(readme_content)

# --- Helpers ---
def safe_string(value):
    if pd.isna(value):
        return None
    value = str(value).strip()
    if not value or value.lower() == "nan":
        return None
    return value

def clean_species_name(value):
    value = safe_string(value)
    return value.lower() if value else None

def clean_class_name(name):
    name = safe_string(name)
    if not name:
        return None
    return re.sub(r"\s*(antibiotic|antibacterial|resistance|agent)?$", "", name.lower()).strip()

def extract_read_species(text):
    if pd.isna(text):
        return []
    entries = str(text).split(";")
    species = []
    for entry in entries:
        clean = re.sub(r"\([^)]*\)", "", entry)
        clean = re.sub(r":.*", "", clean)
        name = safe_string(clean)
        if name:
            species.append(name.lower())
    return species

# --- File Upload for AMR Results ---
amr_file = st.file_uploader("Upload Combined AMR Results CSV", type="csv")

# --- Load MCDA JSON from repo ---
try:
    mcda = pd.read_json("MCDA_matrix_with_alternatives.json")
except Exception as e:
    st.error(f"Error loading MCDA JSON file: {e}")
    st.stop()

if amr_file:
    # --- Load AMR dataset ---
    try:
        df = pd.read_csv(amr_file)
    except Exception as e:
        st.error(f"Error reading uploaded CSV: {e}")
        st.stop()

    # --- Validate AMR file columns ---
    required_amr_cols = {"sample_name", "drug_class", "rpm", "read_species"}
    missing_amr_cols = required_amr_cols - set(df.columns)
    if missing_amr_cols:
        st.error(f"Uploaded AMR Results CSV is missing required columns: {', '.join(sorted(missing_amr_cols))}")
        st.stop()
    else:
        st.success("AMR Results CSV validation passed!")
        st.dataframe(df.head(10))

    # --- Validate MCDA columns ---
    required_mcda_cols = {
        "species", "drug_class", "mortality_score", "incidence_score", "non_fatal_burden_score",
        "transmissibility_score", "preventability_score", "treatability_score",
        "resistance_trend_score", "pipeline_score"
    }
    missing_mcda_cols = required_mcda_cols - set(mcda.columns)
    if missing_mcda_cols:
        st.error(f"MCDA JSON is missing required columns: {', '.join(sorted(missing_mcda_cols))}")
        st.stop()
    else:
        st.success("MCDA JSON validation passed!")
        st.dataframe(mcda.head(10))

    # --- Preprocessing ---
    mcda["species_clean"] = mcda["species"].apply(clean_species_name)
    mcda["drug_class_clean"] = mcda["drug_class"].apply(clean_class_name)

    alt_class_cols = [col for col in mcda.columns if col.startswith("alternative_drug_class")]
    for col in alt_class_cols:
        mcda[col] = mcda[col].apply(clean_class_name)

    score_columns = [
        col for col in mcda.columns
        if col not in ['species', 'drug_class', 'species_clean', 'drug_class_clean', 'all_drug_classes'] + alt_class_cols
    ]

    # Keep numeric score columns only
    score_columns = [col for col in score_columns if pd.api.types.is_numeric_dtype(mcda[col])]

    df["dataset"] = "SingleDataset"
    df["sample_name"] = df["sample_name"].apply(safe_string)
    df["drug_class"] = df["drug_class"].apply(safe_string)
    df["rpm"] = pd.to_numeric(df["rpm"], errors="coerce").fillna(0)

    read_species_extracted = df["read_species"].apply(extract_read_species)
    max_reads = read_species_extracted.apply(len).max()

    if pd.isna(max_reads) or max_reads == 0:
        st.warning("No valid species could be extracted from the read_species column.")
        st.stop()

    max_reads = int(max_reads)

    for i in range(max_reads):
        df[f"read_species_{i+1:02d}"] = read_species_extracted.apply(lambda x: x[i] if i < len(x) else None)

    read_species_cols = [col for col in df.columns if col.startswith("read_species_")]

    for col in read_species_cols:
        df[col] = df[col].apply(clean_species_name)

    # --- Risk Scoring ---
    drug_class_data = []
    detailed_mcda_data = []

    for _, row in df.iterrows():
        rpm = row.get("rpm", 0)
        if pd.isna(rpm) or rpm < 0:
            continue

        if all(pd.isna(row[col]) or row[col] is None for col in read_species_cols):
            continue

        sample = safe_string(row.get("sample_name"))
        if not sample:
            continue

        amr_drug_class = clean_class_name(row.get("drug_class"))
        if not amr_drug_class:
            continue

        for col in read_species_cols:
            species_clean = clean_species_name(row.get(col))
            if not species_clean:
                continue

            mcda_rows_for_species = mcda[mcda["species_clean"] == species_clean]
            if mcda_rows_for_species.empty:
                continue

            direct_match = mcda_rows_for_species[mcda_rows_for_species["drug_class_clean"] == amr_drug_class]
            matched = False
            matched_row = None

            if not direct_match.empty:
                matched_row = direct_match.iloc[0]
                matched = True
            else:
                for _, mcda_row in mcda_rows_for_species.iterrows():
                    alt_classes = [
                        mcda_row[c] for c in alt_class_cols
                        if c in mcda_row and pd.notna(mcda_row[c]) and mcda_row[c] is not None
                    ]
                    if amr_drug_class in alt_classes:
                        matched_row = mcda_row
                        matched = True
                        break

            if not matched or matched_row is None:
                continue

            score_dict = {}
            for score_col in score_columns:
                score_val = pd.to_numeric(matched_row.get(score_col), errors="coerce")
                score_dict[score_col] = 0 if pd.isna(score_val) else score_val

            total_score = sum(score_dict.values()) * rpm

            drug_class_data.append({
                "sample_name": sample,
                "drug_class": amr_drug_class,
                "risk_score": total_score
            })

            detailed_entry = {
                "sample_name": sample,
                "species": species_clean,
                "drug_class": amr_drug_class,
                "rpm": rpm
            }
            detailed_entry.update(score_dict)
            detailed_mcda_data.append(detailed_entry)

    if not drug_class_data or not detailed_mcda_data:
        st.warning("No matched AMR-MCDA records were found after cleaning and filtering.")
        st.stop()

    drug_df = pd.DataFrame(drug_class_data)
    detailed_mcda_df = pd.DataFrame(detailed_mcda_data)

    if drug_df.empty or detailed_mcda_df.empty:
        st.warning("No risk scores could be computed.")
        st.stop()

    drug_class_scores = drug_df.groupby(["sample_name", "drug_class"], as_index=False)["risk_score"].mean()

    all_drug_classes = sorted([x for x in mcda["drug_class_clean"].dropna().unique() if x is not None])

    heatmap_df = drug_class_scores.pivot(index="sample_name", columns="drug_class", values="risk_score").fillna(0)

    for dc in all_drug_classes:
        if dc not in heatmap_df.columns:
            heatmap_df[dc] = 0

    ordered_cols = sorted(heatmap_df.columns.tolist())
    heatmap_df = heatmap_df[ordered_cols]

    heatmap_df["Total_Risk_Score"] = heatmap_df.sum(axis=1)
    heatmap_df["Average_Risk_Score"] = heatmap_df[ordered_cols].mean(axis=1)
    log_heatmap_df = heatmap_df[ordered_cols].apply(lambda col: np.log10(col + 1))

    st.subheader("Log-Transformed AMR Risk Heatmap")
    fig, ax = plt.subplots(figsize=(18, 22))
    sns.heatmap(log_heatmap_df, cmap="YlGnBu", linewidths=0.5, linecolor="gray", ax=ax)
    ax.set_title("Drug Class-Wise Log AMR Risk Scores per Sample", fontsize=16)
    st.pyplot(fig)

    st.subheader("Sample-wise AMR Total Risk Score")
    bar_fig, bar_ax = plt.subplots(figsize=(12, 6))
    heatmap_df["Total_Risk_Score"].sort_values(ascending=False).plot(kind="bar", ax=bar_ax)
    bar_ax.set_ylabel("Total Risk Score")
    bar_ax.set_xlabel("Sample Name")
    bar_ax.set_title("Total AMR Risk Score per Sample")
    bar_ax.tick_params(axis="x", rotation=90)
    st.pyplot(bar_fig)

    # --- Aggregation for Species x Sample x Drug Class ---
    aggregated = detailed_mcda_df.groupby(["sample_name", "species", "drug_class"], as_index=False)["rpm"].sum()
    mcda_scores = detailed_mcda_df.groupby(["species", "drug_class"], as_index=False)[score_columns].mean()
    mcda_scores["mcda_score"] = mcda_scores[score_columns].sum(axis=1)

    final_df = aggregated.merge(
        mcda_scores[["species", "drug_class", "mcda_score"]],
        on=["species", "drug_class"],
        how="left"
    )
    final_df["cumulative_risk_score"] = final_df["rpm"] * final_df["mcda_score"]

    # --- Faceted Heatmaps ---
    st.subheader("Faceted Heatmaps: Species × Sample per Drug Class")
    facet_df = final_df[final_df["species"].notna() & final_df["drug_class"].notna()].copy()

    if not facet_df.empty:
        plot_metric = st.selectbox("Select Metric for Heatmap", ["rpm", "cumulative_risk_score"])
        drug_classes = sorted(facet_df["drug_class"].dropna().unique())

        if drug_classes:
            n_cols = 3
            n_rows = -(-len(drug_classes) // n_cols)

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows), squeeze=False)

            for i, drug_class in enumerate(drug_classes):
                row_i, col_i = divmod(i, n_cols)
                ax = axes[row_i][col_i]
                sub_df = facet_df[facet_df["drug_class"] == drug_class]
                pivot = sub_df.pivot_table(index="species", columns="sample_name", values=plot_metric, aggfunc="sum").fillna(0)

                if not pivot.empty:
                    sns.heatmap(pivot, ax=ax, cmap="YlGnBu", linewidths=0.3, linecolor="gray")
                    ax.set_title(f"Drug Class: {drug_class}", fontsize=12)
                    ax.set_xlabel("Sample")
                    ax.set_ylabel("Species")
                    ax.tick_params(axis="x", labelrotation=90)
                else:
                    ax.axis("off")

            for i in range(len(drug_classes), n_rows * n_cols):
                row_i, col_i = divmod(i, n_cols)
                axes[row_i][col_i].axis("off")

            plt.tight_layout()
            fig.suptitle("Faceted Heatmaps: Species × Sample per Drug Class", fontsize=16, y=1.02)
            st.pyplot(fig)
    else:
        st.info("No faceted heatmap data available after filtering.")

    # --- Trinity-style Aggregation for Export ---
    detailed_agg = (
        detailed_mcda_df
        .groupby(["sample_name", "species", "drug_class"], as_index=False)
        .agg({**{col: "mean" for col in score_columns}, "rpm": "sum"})
    )
    detailed_agg["mcda_score"] = detailed_agg[score_columns].sum(axis=1)
    detailed_agg["cumulative_risk_score"] = detailed_agg["rpm"] * detailed_agg["mcda_score"]

    # --- Downloads ---
    st.subheader("Download Results")
    csv_detailed_mcda = detailed_agg.to_csv(index=False).encode("utf-8")
    csv_heatmap = heatmap_df.reset_index().to_csv(index=False).encode("utf-8")
    csv_final_df = final_df.to_csv(index=False).encode("utf-8")

    st.download_button(
        "Download Detailed MCDA Data as CSV",
        csv_detailed_mcda,
        "Detailed_MCDA_Data.csv",
        "text/csv"
    )
    st.download_button(
        "Download Risk Scores Pivot Data as CSV",
        csv_heatmap,
        "Risk_Scores_Pivot.csv",
        "text/csv"
    )
    st.download_button(
        "Download Aggregated Final Risk Table as CSV",
        csv_final_df,
        "Final_Aggregated_Risk_Scores.csv",
        "text/csv"
    )
