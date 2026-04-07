# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

iDEP (Integrated Differential Expression & Pathway analysis) is a bioinformatics platform for analyzing gene expression data from RNA-Seq, microarray, proteomics, and other methods. It's built as an R package using the Golem framework for Shiny applications.

## Key Commands

### Development
- **Start development server**: Run `source("dev/run_dev.R")` from R console (re-documents and reloads package automatically)
- **Run app in production mode**: `idepGolem::run_app()` from R console
- **Run tests**: `testthat::test_dir("tests/")` or `devtools::test()` or for one specific test only `devtools::test(filter = "deg_full_rank")`
- **Document and reload**: `golem::document_and_reload()` (used in dev workflow)
- **Lint**: uses `.lintr` config enforcing `snake_case` object names

### Package Management
- **Install from GitHub**: `devtools::install_github("https://github.com/gexijin/idepGolem", upgrade = "never")`
- **Install dependencies**: Run the above command first (355 dependencies)

### Docker (Alternative deployment)
- **Run with Docker**: `docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest`
- **Access**: `http://localhost:3838`

## Architecture

### Golem Framework Structure
- **Main entry point**: `R/run_app.R` — loads `idep_data` database once at startup and passes it to every module
- **UI**: `R/app_ui.R` — defines the navbar-based interface structure
- **Server**: `R/app_server.R` — coordinates all modules, manages tab visibility, and passes reactive data downstream
- **Configuration**: `inst/golem-config.yml` manages environment settings

### Shiny Module System
The app has 14 analysis modules (numbered by tab order):
1. `mod_01_load_data` — Data loading and validation
2. `mod_02_pre_process` — Data preprocessing and filtering
3. `mod_03_clustering` — Hierarchical clustering analysis
4. `mod_04_pca` — Principal component analysis
5. `mod_05_deg` — Differential expression analysis (has 2 UI components)
6. `mod_06_pathway` — Pathway enrichment analysis
7. `mod_07_genome` — Genome-wide analysis
8. `mod_08_bicluster` — Biclustering analysis
9. `mod_09_network` — Gene network analysis
10. `mod_10_doc` — Documentation and help
11. `mod_11_enrichment` — Gene set enrichment
12. `mod_12_heatmap` — Heatmap visualization
13. `mod_13_gene_plot` — Individual gene expression plots
14. `mod_14_survey` — User survey

Each module follows the pattern: `mod_XX_name_ui()` and `mod_XX_name_server()`

### Data Flow
The pipeline is hierarchical and one-directional:

```
load_data → pre_process → [clustering, pca, deg, pathway, genome, bicluster, network, heatmap, ...]
```

- Each module **returns a named list of reactive expressions** (not just side effects). Downstream modules receive this list and call its elements as functions.
- The global `idep_data` object (species/pathway databases) is loaded once in `run_app()` and passed as an argument to every module — not reloaded per session.
- A `tab <- reactive(input$navbar)` is passed to modules so they can defer expensive computations until the user visits their tab.
- `app_server.R` contains two large `observe()` blocks controlling tab visibility based on data file format (1=counts, 2=normalized, 3=fold-change+p, 4=fold-change only) and species selection.

### fct_ vs mod_ Separation
- `fct_XX_*.R` files contain **pure functions** with no Shiny dependencies — they take data as arguments and return data. These are unit-testable without a reactive context.
- `mod_XX_*.R` files contain **UI + reactive orchestration** — they call `fct_` functions inside `reactive()` or `eventReactive()` blocks.
- Tests should call `fct_` functions directly (see `tests/testthat/` for examples).

### Key Files
- **Data processing**: `R/fct_01_load_data.R` through `R/fct_13_gene_plot.R`
- **Database utilities**: `R/fct_database.R`
- **UI utilities**: `R/utils_ui.R`
- **Color palettes**: `R/aaa_palette_utils.R` (prefixed `aaa_` to ensure early loading)
- **KEGG pathway rendering**: `R/utils_kegg_pathview.R` — internal replacements for the pathview Bioconductor package (removed as a dependency). Entry point is `mypathview()`, called from `fct_06_pathway.R`. Species validation uses `KEGGREST::keggInfo()` instead of pathview's static `korg` data.

### Global Session-Shared State
`run_app.R` initializes these globals once via Shiny's `onStart` (shared across all user sessions):
- `db_ver` — database version string (e.g. `"data113"`)
- `DATAPATH` — resolved path to the database directory
- `org_info_file` — path to `orgInfo.db`
- `idep_data` — species/demo data loaded by `get_idep_data()`

### Input Data Format Types
The app handles 4 input data types (set in `load_data` module), which drive tab visibility in `app_server.R`:
1. Raw read counts — all tabs visible
2. Normalized expression — all tabs visible
3. Fold-change + FDR — hides Prep, Cluster, PCA, Bicluster, Network for single-comparison; hides PCA, Bicluster, Network for multi-comparison
4. Fold-change only (no FDR) — same as type 3 but also hides Volcano Plot

### Downloadable Reports
RMarkdown templates in `inst/app/www/RMD/` generate downloadable HTML workflow reports (pre-process, clustering, PCA, DEG, pathway). The Prep tab bundles parameters into an `.RData` file so reports can be regenerated locally.

### Database Integration
- Uses environment variable `IDEP_DATABASE` or falls back to `../../data` relative path, then `./data113`
- Database version controlled via `db_ver` variable (currently "data113")
- Species and pathway data loaded via `get_idep_data()` function

### Testing Patterns
- Test data lives in `tests/testthat/testdata/` (synthetic count matrices and design files)
- Tests reference GitHub issues in comments to document the bug being tested (e.g., issue #831, #856)
- `testthat` edition 3 (specified in DESCRIPTION)

### Deployment
- GitHub Actions workflows in `.github/workflows/` build Electron desktop apps for Windows, Mac, and Linux
- Electron builds bundle a full R runtime + all dependencies (~1 hour build time)
- Workflows are manually triggered (workflow_dispatch)

### Dependencies
The app requires 355 R packages including major bioinformatics libraries:
- Bioconductor packages: DESeq2, limma, GSVA, ReactomePA, KEGGREST, KEGGgraph
- Shiny ecosystem: shiny, plotly, DT, visNetwork
- Analysis packages: WGCNA, PCAtools, fgsea, gage
