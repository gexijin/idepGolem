# testdata

Fixture files for `testthat` unit tests and manual UI testing of iDEP.

---

## Expression matrices

### `counts_collinear.csv`

**Format:** Raw integer read counts (format = 1)
**Dimensions:** 200 genes × 9 samples
**Samples:** `A_batch1_1/2/3`, `B_batch2_1/2/3`, `C_batch3_1/2/3`
**Design:** 3 conditions (A, B, C) × 3 replicates, with a confounded batch structure
(condition and batch are perfectly collinear: A↔batch1, B↔batch2, C↔batch3).
**Used by:** Tests that exercise `find_contrast_samples()` with interaction terms
(`I:condition_batch`), `deg_limma()` with a block factor, and any test that
needs a counts matrix with a multi-factor design file.

### `normalized_expr.csv`

**Format:** Log2-scale normalized expression (format = 2, e.g. FPKM/TPM)
**Dimensions:** 200 genes × 9 samples
**Samples:** `A_1, A_2, A_3, B_1, B_2, B_3, C_1, C_2, C_3`
**Design:** 3 conditions (A, B, C) × 3 replicates, single-factor, no batch confounding.
Values are `log2(counts_collinear + 1)` with sample columns renamed so that
`detect_groups()` resolves A_1→A, B_1→B, etc.
**Used by:** All Section A–F tests in `test-normalized_workflow.R` that verify
format = 2 code paths (sample subsetting in heatmap/volcano, `deg_limma` structure,
`deg_information`, `build_deg_gene_table`, etc.).
**Note:** Gene IDs (`gene0001`…`gene0200`) are synthetic and will not match any
species database. Use "Do not convert gene IDs" or a custom GMT file for
end-to-end pathway testing.

---

## Design / sample-information files

The app and the R unit-test layer expect **opposite orientations** for design files.
Two files cover the `normalized_expr.csv` samples; one covers `counts_collinear.csv`.

### `sample_info_simple.csv` — unit-test orientation

**Orientation:** rows = samples, columns = factors
**Samples:** `A_1 … C_3` (9 rows), one column `condition` (values A / B / C)
**Read with:** `read.csv(path, row.names = 1)` → 9 × 1 data frame
**Used by:** R unit tests that call `find_contrast_samples()`,
`deg_limma()`, etc. directly. These functions receive a data frame whose
row names are sample names and whose columns are factor names.
**Do NOT upload to the app** — the app parser expects the transposed format (below).

### `sample_info_ui.csv` — app / UI orientation

**Orientation:** rows = factors, columns = samples
**Factors:** `condition` (one row), values A / B / C across the 9 sample columns
**Read with:** `read.csv(path, row.names = 1)` → 1 × 9 data frame
**Used by:** Manual UI testing. Upload this file as the "Experiment Design" file
in the Load Data tab when using `normalized_expr.csv`.
**Why separate?** The app's design-file parser validates that `colnames(design)`
match `colnames(expression)`, requiring samples as columns. Unit-test helper
functions use the transposed form (`sample_info_simple.csv`) directly.

### `sample_info_collinear.csv` — app / UI orientation, multi-factor

**Orientation:** rows = factors, columns = samples
**Factors (rows):** `condition`, `batch`, `condition_batch` (interaction label)
**Samples:** `A_batch1_1 … C_batch3_3` (9 columns, matching `counts_collinear.csv`)
**Used by:** Manual UI testing with `counts_collinear.csv` (format = 1, counts).
Also referenced by unit tests that pass a sample_info object with `condition` and
`batch` columns (constructed inline, not read from this file directly).

---

## Quick reference

| File | Format | Orientation | Expression file | Usage |
|------|--------|-------------|-----------------|-------|
| `counts_collinear.csv` | counts (1) | — | — | unit tests + UI |
| `normalized_expr.csv` | normalized (2) | — | — | unit tests + UI |
| `sample_info_simple.csv` | design | rows = samples | `normalized_expr.csv` | unit tests only |
| `sample_info_ui.csv` | design | rows = factors | `normalized_expr.csv` | UI manual testing |
| `sample_info_collinear.csv` | design | rows = factors | `counts_collinear.csv` | UI manual testing |
