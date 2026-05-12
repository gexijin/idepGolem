# Manual Testing Checklist: Normalized Expression Data (format = 2)

> **33 items** — Run in order: Load Data → Pre-Process → DEG Stats → DEG (DEG2) → Pathway → Edge Cases.

#### ** All pass as of 3/24/26 in branch v2.4.4_deployed
---

## Test Data Files
#### All files are in tests/testthat/testdata/normalized_expr_testing
| Test Section | Expression File | Design File |
|---|---|---|
| Tab 1–6, E2, E3 | `normalized_expr.csv` | `sample_info_normalized.csv` |
| Tab 7 (7.1–7.4) | iDEP built-in demo data (see note below) | loaded automatically with demo |
| 7.5 (real gene symbols) | `real_gene_sample_info_normalized.csv` | `sample_info_normalized.csv` |
| E1 | Start: `normalized_expr.csv`, then iDEP demo count data | `sample_info_normalized.csv` first |
| E4 (2-sample) | `normalized_expr_2_sample.csv` | `sample_info_normalized_2_samp.csv` |

All files are in `tests/testthat/testdata/`.

### Design file format — important
Always upload **`sample_info_normalized.csv`** — it uses the format iDEP expects: **samples as columns, factors as rows** (the "transposed" layout). Do NOT use any design file where samples are rows; the app will reject it with:
> *"Sample information file not recognized. Column names must be exactly the same. Each row is a factor. Each column represent a sample."*

### Expression file note
`normalized_expr.csv` contains **50 genes** with built-in differential expression:
- **Genes 1–10**: UP in B vs A and B vs C
- **Genes 11–20**: DOWN in B vs A and B vs C
- **Genes 21–30**: UP in C vs A and C vs B
- **Genes 31–40**: DOWN in C vs A and C vs B
- **Genes 41–50**: background (no change across groups)

This structure makes it possible to test all DEG tab subtabs without real biological data.

### Demo data note (Tab 7)
On the Load Data tab, click **"Load demo data"** (or equivalent button) to load iDEP's built-in example dataset. This dataset uses real Ensembl IDs that map to KEGG pathways. The fake IDs in `normalized_expr.csv` (`gene0001`, `gene0002`, etc.) will return zero pathway hits and should **not** be used for Tab 7 tests 7.1–7.4.

---

## Tab 1: Load Data

> Files: `normalized_expr.csv` + `sample_info_normalized.csv`

- [ ] **1.1** Select **"Normalized expression"** from the data type dropdown and upload `normalized_expr.csv` — verify the file is accepted without error

- [ ] **1.3** Upload `sample_info_normalized.csv` as the design file — verify the experiment design table populates with 9 samples (A_1–A_3, B_1–B_3, C_1–C_3) assigned to conditions A, B, C

- [ ] **1.4** Check **"Do not convert gene IDs"** — verify it is selectable and the selection persists after navigating away and back

---

## Tab 2: Pre-Process

> Files: `normalized_expr.csv` + `sample_info_normalized.csv`. Continue from Tab 1.

- [ ] **2.4** Set min expression level and min samples — verify the gene count updates in response

- [ ] **2.5** Click **"Processed Data"** download — verify the CSV opens correctly and contains sample columns A_1 through C_3

- [ ] **2.6** Click **Pre-Process report** download — verify the HTML report generates without error

---

## Tab 5: DEG Stats (DEG1)

> Files: `normalized_expr.csv` + `sample_info_normalized.csv`
>
> ⚠️ **Test 5.4 must be done before the design file is loaded.** If you already uploaded the design file in Tab 1, either start a fresh session for 5.4 or clear the design file first, then re-upload it before 5.5.

- [ ] **5.4** With **no design file loaded**: verify groups are auto-detected from column names (A_1/A_2/A_3 → group A, B_1/B_2/B_3 → group B, etc.) and Submit runs without error

- [ ] **5.5** With **`sample_info_normalized.csv` loaded**: verify the factor dropdown shows "condition", the comparison selector populates (expect comparisons such as B-A, C-A, C-B), and Submit runs without error

- [ ] **5.6** After Submit: verify the **Results** subtab renders — expect a bar chart showing DEG counts and a summary table with Up/Down counts per comparison

- [ ] **5.7** Click **"Results" download** (LFC table) — verify the CSV contains logFC, adjPval, and expression columns

- [ ] **5.8** Click **"Gene Lists" download** — verify the CSV contains only Up and Down rows (None rows should be filtered out)

- [ ] **5.10** Select **2+ comparisons** (e.g., B-A and C-A) and Submit — verify the Venn diagram renders

- [ ] **5.11** Click **DEG Stats report** download — verify the HTML generates and contains "limma" (not "DESeq2") in the method description

---

## Tab 6: DEG (DEG2)

> Files: `normalized_expr.csv` + `sample_info_normalized.csv`
> Run a **B vs A** comparison in DEG Stats first, then navigate to this tab.

### Genes subtab

- [ ] **6.2** Toggle the direction filter between "Up-regulated" and "Down-regulated" — verify the gene table updates each time

- [ ] **6.3** Verify the **"Gene List" download** button is visible in the sidebar when on the Genes subtab

- [ ] **6.6** Click a gene row — verify the expression popup appears

- [ ] **6.7** Verify **no "Raw data" option** appears in the popup dropdown (format=2)

### Heatmap subtab

- [ ] **6.9** Verify the gene number selector and sort method both update the heatmap

### Volcano Plot subtab

- [ ] **6.11** Verify the volcano plot renders with correct up/down coloring

- [ ] **6.12** Verify labeling and download works

### MA Plot subtab

- [ ] **6.13** Verify the MA plot renders

- [ ] **6.14** Verify contrast-only averaging behavior (as described)

### Scatter Plot subtab

- [ ] **6.15** Verify the scatter plot renders and download works

### Report

- [ ] **6.16** Verify DEG report generates with limma label

---

## Tab 7: Pathway

- [ ] **7.1** Comparison dropdown populates

- [ ] **7.2** KEGG results render

- [ ] **7.3** No-ID conversion produces graceful empty result (no crash)

- [ ] **7.4** Downloads open correctly

- [ ] **7.5** Real gene symbols produce valid enrichment results (use B-A comparison in Stats/DEG/Pathway tabs and GSEA pathway method)

---

## Edge Cases

- [ ] **E1** Format switch resets UI and restores DESeq2

- [ ] **E2** No-ID conversion preserves original gene IDs in outputs

- [ ] **E3** Missing comparison produces warning, no crash

- [ ] **E4** Two-sample dataset handled gracefully

---

## Additional Tests (completeness rather than necessary tests)


### Load Data

- [ ] **M1** Upload a normalized file containing **negative values** — verify the file is accepted (no rejection; counts-only validation must NOT trigger)

### Pre-Process (UI Behavior)

- [ ] **M2** Verify the **"Reads" subtab is NOT present** for normalized data

- [ ] **M3** Verify filter labels reflect **normalized units (FPKM/TPM style)**, not CPM

- [ ] **M4** Toggle **log transformation** — verify the distribution plot visibly changes

### DEG Stats (UI Gating)

- [ ] **M5** Verify **DESeq2 is NOT available** as a method option

- [ ] **M6** Verify DESeq2-specific options are hidden:
  - "Threshold Wald Test"
  - "Independent filtering"

- [ ] **M7** Verify **reference level dropdowns are NOT shown** when selecting factors

- [ ] **M8** Open the **R Code tab** — verify it contains **limma code only** (no DESeq2)

### DEG (Genes + Downloads)

- [ ] **M9** Verify gene table columns:
  - ID, Ensembl ID, log2 FC, Adj. Pval, Description

- [ ] **M10** Click **Gene List download** — verify CSV contains:
  - `ensembl_ID`, `symbol`, `entrezgene_id`, `log2FC`, `Adjusted_P_value`, `description`

- [ ] **M11** Verify downloaded CSV reflects **current direction filter only** (Up or Down)

### DEG (Heatmap + Data Integrity)

- [ ] **M12** Verify heatmap shows **ONLY contrast samples** (e.g., 6 columns for A vs B, not 9)

- [ ] **M13** Download **Heatmap Data** — verify column count matches heatmap

### MA Plot (UI Validation)

- [ ] **M14** Visually confirm averages are NOT influenced by non-contrast groups (complements 6.14 with a second check)

### Pathway (Additional Coverage)

- [ ] **M15** Verify **gene symbols (not Ensembl IDs)** appear in enrichment results when using real symbol-based input with no conversion

### Edge Case Enhancements

- [ ] **M16** After format switch (E1), verify:
  - No stale plots remain
  - No stale downloads persist

- [ ] **M17** Attempt invalid workflow transitions (e.g., skipping DEG Stats → going to DEG tab) — verify app does not crash and shows appropriate empty state

