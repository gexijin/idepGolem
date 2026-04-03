# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

iDEP (Integrated Differential Expression & Pathway analysis) is a bioinformatics platform for analyzing gene expression data from RNA-Seq, microarray, proteomics, and other methods. It's built as an R package using the Golem framework for Shiny applications.

## Key Commands

### Development
- **Start development server**: Run `source("dev/run_dev.R")` from R console
- **Run app in production mode**: `idepGolem::run_app()` from R console
- **Run tests**: `testthat::test_dir("tests/")` or `devtools::test()` or for one specific test only `devtools::test(filter = "deg_full_rank")`
- **Document and reload**: `golem::document_and_reload()` (used in dev workflow)

### Package Management
- **Install from GitHub**: `devtools::install_github("https://github.com/gexijin/idepGolem", upgrade = "never")`
- **Install dependencies**: Run the above command first (355 dependencies)

### Docker (Alternative deployment)
- **Run with Docker**: `docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest`
- **Access**: `http://localhost:3838`

## Architecture

### Golem Framework Structure
- **Main entry point**: `R/run_app.R` exports the main `run_app()` function
- **UI**: `R/app_ui.R` defines the navbar-based interface structure
- **Server**: `R/app_server.R` coordinates all modules and global settings
- **Configuration**: `inst/golem-config.yml` manages environment settings

### Shiny Module System
The app is organized into 12 main analysis modules:
1. `mod_01_load_data` - Data loading and validation
2. `mod_02_pre_process` - Data preprocessing and filtering
3. `mod_03_clustering` - Hierarchical clustering analysis
4. `mod_04_pca` - Principal component analysis
5. `mod_05_deg` - Differential expression analysis (has 2 UI components)
6. `mod_06_pathway` - Pathway enrichment analysis
7. `mod_07_genome` - Genome-wide analysis
8. `mod_08_bicluster` - Biclustering analysis
9. `mod_09_network` - Gene network analysis
10. `mod_10_doc` - Documentation and help
11. `mod_11_enrichment` - Gene set enrichment
12. `mod_12_heatmap` - Heatmap visualization

Each module follows the pattern: `mod_XX_name_ui()` and `mod_XX_name_server()`

### Data Flow
- Data passes between modules using reactive expressions
- `load_data` module output feeds into `pre_process`
- Processed data flows to downstream analysis modules
- Global configuration in `app_server.R` sets database paths and options

### Key Functions by Module
- **Data processing**: `R/fct_01_load_data.R` through `R/fct_12_heatmap.R`, plus `R/fct_13_gene_plot.R`
- **Database utilities**: `R/fct_database.R`
- **UI utilities**: `R/utils_ui.R`
- **Analysis utilities**: `R/utils_analysis_random.R`

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

### Dependencies
The app requires 355 R packages including major bioinformatics libraries:
- Bioconductor packages: DESeq2, limma, GSVA, ReactomePA
- Shiny ecosystem: shiny, plotly, DT, visNetwork
- Analysis packages: WGCNA, PCAtools, fgsea, gage