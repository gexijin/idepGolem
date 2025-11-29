#!/usr/bin/env bash
set -euo pipefail
VER=4.4.1
OUT="r-linux"
mkdir -p "$OUT"
cd "$OUT"

# Prebuilt Linux tarball (x86_64)
# If needed, switch to distribution-specific builds
curl -LO https://cloud.r-project.org/src/base/R-4/R-${VER}.tar.gz
tar -xzf R-${VER}.tar.gz
# Build minimal R (requires build deps on the runner) – or swap to prebuilt if you prefer
cd R-${VER}
./configure --prefix="$(pwd)/../R"
make -j2
make install
cd ..
rm -rf "R-${VER}" "R-${VER}.tar.gz"
echo "Linux R installed under r-linux/R"







----------- DEV ONLY testing

a) fetch DB
#!/usr/bin/env Rscript
# fetch_db.R — optional helper for dev/admin use

args <- commandArgs(trailingOnly = TRUE)
db_ver <- Sys.getenv("SHINYGO_DB_VER", if (length(args) >= 1) args[[1]] else "data113")
base   <- Sys.getenv("SHINYGO_DB_URL", if (length(args) >= 2) args[[2]] else "http://bioinformatics.sdstate.edu/data")
parent <- Sys.getenv("IDEP_DATABASE", if (length(args) >= 3) args[[3]] else file.path(Sys.getenv("HOME"), ".shinygo-data"))

dir.create(parent, recursive = TRUE, showWarnings = FALSE)
dest_ver <- file.path(parent, db_ver)
tgz <- file.path(parent, paste0(db_ver, ".tar.gz"))

url <- sprintf("%s/%s/%s.tar.gz", base, db_ver, db_ver)
message("Downloading: ", url)
utils::download.file(url, tgz, mode = "wb", quiet = FALSE, timeout = 300)

message("Extracting to: ", parent)
utils::untar(tgz, exdir = parent)
unlink(tgz)

# sanity check + current marker
org_info <- file.path(dest_ver, "demo", "orgInfo.db")
if (!file.exists(org_info)) {
  stop("orgInfo.db not found at: ", org_info)
}
writeLines(db_ver, file.path(parent, "current"))
message("Done. IDEP_DATABASE parent is: ", parent, "   version: ", db_ver)

###how to use
Rscript fetch_db.R

# explicit version + custom parent folder
IDEP_DATABASE="D:/idep-data" Rscript fetch_db.R data113 "http://bioinformatics.sdstate.edu/data"