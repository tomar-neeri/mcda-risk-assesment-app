# 🧬 MCDA-AMR: A Multi-Criteria Decision Analysis Framework for Antimicrobial Resistance Risk Assessment from Metagenomic Datasets

## 📌 Overview

**MCDA-AMR** is a user-friendly, web-based analytical tool designed to assess and prioritize the risk of **antimicrobial resistance (AMR)** based on **metagenomic data**. It uses a structured **Multi-Criteria Decision Analysis (MCDA)** approach, integrating gene abundance, species-level risk, and drug class severity to provide actionable insights into AMR surveillance and public health prioritization.

---

## 🔬 Scientific Background

Antimicrobial Resistance (AMR) is a critical global public health issue. According to WHO, AMR directly caused 1.27 million deaths in 2019 and contributed to nearly 5 million deaths. Metagenomic sequencing has enabled detection of AMR genes in diverse environments, but interpreting the resulting high-dimensional data in a clinically meaningful way remains a challenge.

### Why MCDA?

Multi-Criteria Decision Analysis (MCDA) allows integration of multiple biological, clinical, and ecological factors—such as **mortality, prevalence, resistance trends, and treatability**—to assess risk. This framework fills the gap between raw sequencing data and public health decision-making.

---

## 🛠️ Features

- Upload AMR gene detection results and MCDA risk matrix
- Compute **species-drug abundance-adjusted risk scores**
- Perform **Trinity aggregation** (Species × Drug × Sample)
- Visualize risk with:
  - Sample-wise bar plots
  - Species-wise heatmaps
  - Log-transformed risk matrices
- Export risk-scored tables for further analysis

---

## 📁 Input Files

### 1. AMR Results File (CSV)

| sample_name | read_species             | drug_class     | rpm   |
|-------------|--------------------------|----------------|-------|
| Sample_001  | Klebsiella pneumoniae    | Carbapenem     | 150.5 |
| Sample_001  | Escherichia coli         | Fluoroquinolone| 230.2 |
| Sample_002  | Acinetobacter baumannii  | Carbapenem     | 85.0  |

---

### 2. MCDA Risk Matrix File (CSV)

| species                | drug_class     | mortality_score | incidence_score | ... | pipeline_score | alternative_drug_class |
|------------------------|----------------|------------------|------------------|-----|----------------|--------------------------|
| Klebsiella pneumoniae  | Carbapenem     | 5                | 4                | ... | 1.5            | Carbapenems              |

Each scoring feature (Mortality, Incidence, etc.) ranges from **0.5 to 5**.

---

## 📥 Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/tomar-neeri/mcda-risk-assesment-app.git 
cd mcda-amr-app
pip install -r requirements.txt
```

If you don’t have Python installed, install it first: https://www.python.org/downloads/

---

## 🚀 Running the App Directly

Launch the Streamlit web app locally:

```bash
streamlit run mcda-amr-app.py
```
OR 

Then open [ https://mcda-amr-risk-assessment-pipeline-eepm-neeri.streamlit.app/ ] in your browser.

---

## 🖼️ Outputs

1. **Triplet Scores Table**  
   Cumulative risk scores for each Sample × Species × Drug combination.

2. **Risk Pivot Table**  
   Summarized risk scores per sample for comparison.

3. **Aggregated Risk Table**  
   Species- and drug class-level cumulative risk metrics.

4. **Visualizations**  
   - Heatmaps of log-transformed risk scores  
   - Sample-wise and species-wise bar plots  
   - RPM vs Risk comparisons

All outputs are downloadable in CSV format from the dashboard.

---

## 🧪 Example Datasets

Example input files are included:
- `sample_AMR_results.csv`
- `sample_MCDA_matrix.csv`

You can use these to test the app or replace them with your own metagenomic data.

---

## 🔁 Risk Score Formula

```text
Risk Score = (Sum of MCDA Attributes) × RPM
```

For triplet aggregation:
```text
Cumulative Risk Score = RPM_sum × Mean_MCDA_Score
```

---

## 📚 References

1. WHO. Global burden of bacterial AMR in 2019: a systematic analysis. Lancet (2022).
2. Wellington et al. The role of the natural environment in the emergence of antibiotic resistance. Environ Int. (2013).
3. Fitzpatrick & Walsh. Antibiotic resistance genes across a wide variety of metagenomes. FEMS Microbiol Rev. (2016).
4. Martínez JL et al. Ranking ARGs based on risk to human health. Nat Rev Microbiol. (2015).
5. Zhang et al. Classifying ARGs using omics and ecological context. ISME J. (2021).
6. Arango-Argoty et al. DeepARG: A deep learning approach for detecting ARGs. Microbiome (2018).
7. Liao et al. PathoFact: Integrated pipeline for AMR, virulence, and toxin detection. Bioinformatics (2021).
8. CARD: Comprehensive Antibiotic Resistance Database. https://card.mcmaster.ca/
9. Chan Zuckerberg ID Platform: https://czid.org/
10. WHO Bacterial Priority Pathogens List 2024.

---

## 🧑‍💻 Authors & Contributions

Developed by Siddharth Singh Tomar  
Data: SARS-CoV-2 genomic surveillance, Central India, March–April 2023  
Framework validated on 96 upper respiratory metagenomes using CZID and CARD v3.2.6

---

## 📜 License

MIT License © 2025 Siddharth Singh Tomar

---

## 💬 Contact

For queries, feedback, or collaborations, feel free to reach out at:  
📧 youremail@example.com  
🔗 https://github.com/tomar-neeri/mcda-risk-assesment-app.git 