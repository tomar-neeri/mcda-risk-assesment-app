import streamlit as st
import pandas as pd
import numpy as np
import re
import seaborn as sns
import matplotlib.pyplot as plt
from io import BytesIO

st.set_page_config(layout="wide")
st.title("A MCDA Framework for Assessing and Prioritizing AMR Risks from Metagenomic Datasets")

# --- File Upload ---
amr_file = st.file_uploader("Upload Combined AMR Results CSV", type="csv")
mcda_file = st.file_uploader("Upload MCDA Matrix CSV", type="csv")

if amr_file and mcda_file:
    # --- Load datasets ---
    df = pd.read_csv(amr_file)
    mcda = pd.read_csv(mcda_file)

    # --- Validate AMR file columns ---
    required_amr_cols = {"sample_name", "drug_class", "rpm", "read_species"}
    missing_amr_cols = required_amr_cols - set(df.columns)
    if missing_amr_cols:
        st.error(f"Uploaded AMR Results CSV is missing required columns: {', '.join(missing_amr_cols)}")
        st.stop()
    else:
        st.success("AMR Results CSV validation passed!")
        st.dataframe(df.head(10))

    # --- Validate MCDA file columns ---
    required_mcda_cols = {
        "species", "drug_class", "mortality_score", "incidence_score", "non_fatal_burden_score",
        "transmissibility_score", "preventability_score", "treatability_score",
        "resistance_trend_score", "pipeline_score"
    }
    missing_mcda_cols = required_mcda_cols - set(mcda.columns)
    if missing_mcda_cols:
        st.error(f"Uploaded MCDA Matrix CSV is missing required columns: {', '.join(missing_mcda_cols)}")
        st.stop()
    else:
        st.success("MCDA Matrix CSV validation passed!")
        st.dataframe(mcda.head(10))

    # --- Preprocessing ---
    mcda["species_clean"] = mcda["species"].str.strip().str.lower()

    def clean_class_name(name):
        if pd.isna(name):
            return None
        return re.sub(r"\s*(antibiotic|antibacterial|resistance|agent)?$", "", name.strip().lower())

    mcda["drug_class_clean"] = mcda["drug_class"].apply(clean_class_name)
    alt_class_cols = [col for col in mcda.columns if col.startswith("alternative_drug_class")]
    for col in alt_class_cols:
        mcda[col] = mcda[col].apply(clean_class_name)

    score_columns = [
        col for col in mcda.columns if col not in ['species', 'drug_class', 'species_clean'] + alt_class_cols + ['all_drug_classes', 'drug_class_clean']
    ]

    def extract_read_species(text):
        if pd.isna(text):
            return []
        entries = text.split(';')
        species = []
        for entry in entries:
            clean = re.sub(r'\([^)]*\)', '', entry)
            clean = re.sub(r':.*', '', clean)
            name = clean.strip()
            if name:
                species.append(name.lower())
        return species

    df["dataset"] = "SingleDataset"
    read_species_extracted = df["read_species"].apply(extract_read_species)
    max_reads = read_species_extracted.apply(len).max()
    for i in range(max_reads):
        df[f'read_species_{i+1:02d}'] = read_species_extracted.apply(lambda x: x[i] if i < len(x) else None)

    read_species_cols = [col for col in df.columns if col.startswith("read_species")]
    for col in read_species_cols:
        df[col] = df[col].astype(str).str.strip().str.lower()

    # --- Risk Scoring ---
    drug_class_data = []
    detailed_mcda_data = []

    for _, row in df.iterrows():
        rpm = row.get("rpm", 0)
        if rpm < 0:
            continue
        if all(pd.isna(row[col]) or str(row[col]).strip().lower() == "nan" for col in read_species_cols):
            continue
        sample = row["sample_name"]
        amr_drug_class_raw = row.get("drug_class")
        if not isinstance(amr_drug_class_raw, str) or not amr_drug_class_raw.strip():
            continue
        amr_drug_class = clean_class_name(amr_drug_class_raw)

        for col in read_species_cols:
            species = row.get(col)
            if not species or species == "nan":
                continue
            species_clean = species.strip().lower()
            mcda_rows_for_species = mcda[mcda["species_clean"] == species_clean]
            direct_match = mcda_rows_for_species[mcda_rows_for_species["drug_class_clean"] == amr_drug_class]
            matched = False
            if not direct_match.empty:
                matched_row = direct_match.iloc[0]
                matched = True
            else:
                for _, mcda_row in mcda_rows_for_species.iterrows():
                    alt_classes = [mcda_row[col] for col in alt_class_cols if pd.notna(mcda_row[col])]
                    if amr_drug_class in alt_classes:
                        matched_row = mcda_row
                        matched = True
                        break
            if not matched:
                continue
            score_dict = {col: matched_row[col] for col in score_columns}
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

    drug_df = pd.DataFrame(drug_class_data)
    detailed_mcda_df = pd.DataFrame(detailed_mcda_data)

    drug_class_scores = drug_df.groupby(["sample_name", "drug_class"])["risk_score"].mean().reset_index()
    all_drug_classes = sorted(set(mcda["drug_class_clean"].dropna().unique()))
    heatmap_df = drug_class_scores.pivot(index="sample_name", columns="drug_class", values="risk_score").fillna(0)

    for dc in all_drug_classes:
        if dc not in heatmap_df.columns:
            heatmap_df[dc] = 0
    ordered_cols = sorted(heatmap_df.columns)
    heatmap_df = heatmap_df[ordered_cols]

    heatmap_df["Total_Risk_Score"] = heatmap_df.sum(axis=1)
    heatmap_df["Average_Risk_Score"] = heatmap_df[ordered_cols].mean(axis=1)
    log_heatmap_df = heatmap_df[ordered_cols].applymap(lambda x: np.log10(x + 1))

    st.subheader("Log-Transformed AMR Risk Heatmap")
    fig, ax = plt.subplots(figsize=(18, 22))
    sns.heatmap(log_heatmap_df, cmap="YlGnBu", linewidths=0.5, linecolor='gray', ax=ax)
    plt.title("Drug Class-Wise Log AMR Risk Scores per Sample", fontsize=16)
    st.pyplot(fig)

    st.subheader("Sample-wise AMR Total Risk Score")
    bar_fig, bar_ax = plt.subplots(figsize=(12, 6))
    heatmap_df["Total_Risk_Score"].sort_values(ascending=False).plot(kind='bar', ax=bar_ax)
    bar_ax.set_ylabel("Total Risk Score")
    bar_ax.set_xlabel("Sample Name")
    bar_ax.set_title("Total AMR Risk Score per Sample")
    bar_ax.tick_params(axis='x', rotation=90)
    st.pyplot(bar_fig)

    # --- Fixed Aggregation for Species x Sample x Drug Class ---
    aggregated = detailed_mcda_df.groupby(['sample_name', 'species', 'drug_class'], as_index=False)['rpm'].sum()

    mcda_scores = detailed_mcda_df.groupby(['species', 'drug_class'])[score_columns].mean().reset_index()
    mcda_scores['mcda_score'] = mcda_scores[score_columns].sum(axis=1)

    final_df = aggregated.merge(mcda_scores[['species', 'drug_class', 'mcda_score']], on=['species', 'drug_class'], how='left')
    final_df['cumulative_risk_score'] = final_df['rpm'] * final_df['mcda_score']

    # --- Faceted Heatmaps ---
    st.subheader("Faceted Heatmaps: Species × Sample per Drug Class")
    facet_df = final_df[final_df['species'].notna() & final_df['drug_class'].notna()].copy()
    plot_metric = st.selectbox("Select Metric for Heatmap", ["rpm", "cumulative_risk_score"])
    drug_classes = sorted(facet_df['drug_class'].unique())
    n_cols = 3
    n_rows = -(-len(drug_classes) // n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows), squeeze=False)

    for i, drug_class in enumerate(drug_classes):
        row, col = divmod(i, n_cols)
        ax = axes[row][col]

        sub_df = facet_df[facet_df["drug_class"] == drug_class]
        pivot = sub_df.pivot_table(index="species", columns="sample_name", values=plot_metric, aggfunc='sum').fillna(0)

        if not pivot.empty:
            sns.heatmap(pivot, ax=ax, cmap="YlGnBu", linewidths=0.3, linecolor='gray')
            ax.set_title(f"Drug Class: {drug_class}", fontsize=12)
            ax.set_xlabel("Sample")
            ax.set_ylabel("Species")
            ax.tick_params(axis='x', labelrotation=90)
        else:
            ax.axis('off')

    for i in range(len(drug_classes), n_rows * n_cols):
        row, col = divmod(i, n_cols)
        axes[row][col].axis('off')

    plt.tight_layout()
    fig.suptitle("Faceted Heatmaps: Species × Sample per Drug Class", fontsize=16, y=1.02)
    st.pyplot(fig)

    # --- Trinity-style Aggregation for Detailed MCDA Export ---
    detailed_agg = (
        detailed_mcda_df
        .groupby(["sample_name", "species", "drug_class"], as_index=False)
        .agg({**{col: "mean" for col in score_columns}, "rpm": "sum"})
    )
    detailed_agg["mcda_score"] = detailed_agg[score_columns].sum(axis=1)
    detailed_agg["cumulative_risk_score"] = detailed_agg["rpm"] * detailed_agg["mcda_score"]

    # --- Downloads ---
    st.subheader("Download Results")

    csv_detailed_mcda = detailed_agg.to_csv(index=False).encode('utf-8')
    csv_heatmap = heatmap_df.reset_index().to_csv(index=False).encode('utf-8')
    csv_final_df = final_df.to_csv(index=False).encode('utf-8')

    st.download_button(
        label="📥 Download Detailed MCDA Data as CSV",
        data=csv_detailed_mcda,
        file_name="Detailed_MCDA_Data.csv",
        mime="text/csv"
    )

    st.download_button(
        label="📥 Download Risk Scores Pivot Data as CSV",
        data=csv_heatmap,
        file_name="Risk_Scores_Pivot.csv",
        mime="text/csv"
    )

    st.download_button(
        label="📥 Download Aggregated Final Risk Table as CSV",
        data=csv_final_df,
        file_name="Final_Aggregated_Risk_Scores.csv",
        mime="text/csv"
    )
