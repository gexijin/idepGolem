#!/usr/bin/env bash
set -euo pipefail

# ==================== Config ====================
VER="${R_VERSION:-4.5.2}"   # override with env R_VERSION=...
# Auto-detect ARCH unless overridden via R_ARCH
if [[ -z "${R_ARCH:-}" ]]; then
  case "$(uname -m)" in
    arm64|aarch64) ARCH="arm64" ;;
    x86_64|amd64)  ARCH="x86_64" ;;
    *)             ARCH="x86_64" ;;
  esac
else
  ARCH="${R_ARCH}"
fi

# Candidate URLs (primary + fallbacks)
CANDIDATES=(
  "https://cran.r-project.org/bin/macosx/big-sur-${ARCH}/base/R-${VER}-${ARCH}.pkg"
  "https://cloud.r-project.org/bin/macosx/big-sur-${ARCH}/base/R-${VER}-${ARCH}.pkg"
  # Intel builds are sometimes published without the -x86_64 suffix
  "https://cran.r-project.org/bin/macosx/big-sur-x86_64/base/R-${VER}.pkg"
  "https://cloud.r-project.org/bin/macosx/big-sur-x86_64/base/R-${VER}.pkg"
)

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
OUTDIR="${SCRIPT_DIR}/r-mac"

RFRAMEWORK_DEST="${OUTDIR}/R.framework"
RDEST="${OUTDIR}/R"   # <--- this is what build-electron.yml expects

TMP="$(mktemp -d)"
cleanup(){ rm -rf "$TMP"; }
trap cleanup EXIT

PKG_PATH=""
for URL in "${CANDIDATES[@]}"; do
  echo "Trying ${URL} ..."
  if curl -fL --retry 3 --connect-timeout 20 -o "${TMP}/R.pkg" "${URL}"; then
    PKG_PATH="${TMP}/R.pkg"
    echo "Downloaded: ${URL}"
    break
  fi
done

if [[ -z "${PKG_PATH}" ]]; then
  echo "ERROR: Unable to download R ${VER} pkg for macOS (tried ${#CANDIDATES[@]} URLs)." >&2
  exit 1
fi

echo "Expanding pkg ..."
pkgutil --expand-full "${PKG_PATH}" "${TMP}/expanded"

# Locate R.framework inside expanded pkg
RFW="$(/usr/bin/find "${TMP}/expanded" -type d -name 'R.framework' -print -quit || true)"
if [[ -z "${RFW}" ]]; then
  echo "ERROR: R.framework not found in expanded package:" >&2
  /usr/bin/find "${TMP}/expanded" -maxdepth 4 -print
  exit 1
fi

echo "Copying R.framework to ${RFRAMEWORK_DEST} and flattening into ${RDEST} ..."
rm -rf "${RFRAMEWORK_DEST}" "${RDEST}"
mkdir -p "${OUTDIR}"

# Keep the full framework (optional but nice to have)
ditto "${RFW}" "${RFRAMEWORK_DEST}"

# Flatten framework layout so we get r-mac/R/bin/R, etc.
ditto "${RFRAMEWORK_DEST}/Resources" "${RDEST}"

echo "Rscript version:"
"${RDEST}/bin/Rscript" --version
echo "✅ macOS R runtime ready at: ${RDEST}"
