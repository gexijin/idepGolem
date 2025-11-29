# electron/scripts/get_r_windows.ps1
$ErrorActionPreference = "Stop"

# -------- Config --------
$Rver = $env:R_VERSION
if (-not $Rver -or $Rver -eq "") { $Rver = "4.4.2" }   # default version

Write-Host "Using R version: $Rver"

# runtime destination (flat layout)
# script is run from electron/scripts so ../runtime/win/R is the target
$repoRoot = Resolve-Path ".." | Select-Object -ExpandProperty Path
$destWin  = Join-Path $repoRoot "runtime\win"
$destR    = Join-Path $destWin  "R"
New-Item -ItemType Directory -Force -Path $destWin | Out-Null

# Temp working dir
$tmp = Join-Path $env:TEMP ("rwin_" + [Guid]::NewGuid())
New-Item -ItemType Directory -Force -Path $tmp | Out-Null

try {
  # -------- Download installer --------
  $exe = Join-Path $tmp "R-$Rver-win.exe"
  $urls = @(
    "https://cloud.r-project.org/bin/windows/base/old/$Rver/R-$Rver-win.exe",
    "https://cran.r-project.org/bin/windows/base/old/$Rver/R-$Rver-win.exe",
    "https://cloud.r-project.org/bin/windows/base/R-$Rver-win.exe",
    "https://cran.r-project.org/bin/windows/base/R-$Rver-win.exe"
  )

  $downloaded = $false
  foreach ($u in $urls) {
    Write-Host "Trying $u ..."
    try {
      Invoke-WebRequest -Uri $u -OutFile $exe -UseBasicParsing
      if ((Get-Item $exe).Length -gt 0) { $downloaded = $true; Write-Host "Downloaded $u"; break }
    } catch { Write-Host "Download failed from $u, trying next mirror..." }
  }
  if (-not $downloaded) { throw "Failed to download R-$Rver Windows installer from all candidates." }

  # -------- Install silently into temp base (installer creates R-x.y.z under here) --------
  $installBase = Join-Path $tmp "R-install"
  New-Item -ItemType Directory -Force -Path $installBase | Out-Null

  Write-Host "Installing R silently into $installBase ..."
  & $exe /VERYSILENT /DIR="$installBase" /NORESTART /SP- /SUPPRESSMSGBOXES | Out-Null

  # -------- Locate the versioned R home (R-x.y.z) --------
  # Prefer by discovering Rscript.exe, then take its parent twice to get R-x.y.z
  $rscript = Get-ChildItem -Recurse -Path $installBase -Filter Rscript.exe -File -ErrorAction SilentlyContinue | Select-Object -First 1
  if (-not $rscript) { throw "Rscript.exe not found under $installBase after install." }

  $rHome = $rscript.Directory.Parent  # ...\R-x.y.z
  if (-not (Test-Path (Join-Path $rHome.FullName "bin\R.exe"))) {
    throw "R.exe not found in $(Join-Path $rHome.FullName 'bin'); unexpected layout."
  }

  Write-Host "Detected R home: $($rHome.FullName)"

  # -------- Normalize to flat layout: ../runtime/win/R/{bin,library,...} --------
  if (Test-Path $destR) {
    Write-Host "Cleaning existing $destR ..."
    Remove-Item -Recurse -Force $destR
  }
  New-Item -ItemType Directory -Force -Path $destR | Out-Null

  # Copy CONTENTS of R-x.y.z into ../runtime/win/R (so we get R/bin, not R/R-x.y.z/bin)
  Write-Host "Copying portable R to $destR ..."
  Copy-Item -Recurse -Force -Path (Join-Path $rHome.FullName "*") -Destination $destR

  # -------- Sanity checks --------
  $destRscript = Join-Path $destR "bin\Rscript.exe"
  if (-not (Test-Path $destRscript)) { throw "Missing $destRscript after copy." }
  $destRexe = Join-Path $destR "bin\R.exe"
  if (-not (Test-Path $destRexe)) { throw "Missing $destRexe after copy." }
  $destRdll = Join-Path $destR "bin\x64\R.dll"
  if (-not (Test-Path $destRdll)) { Write-Host "Warning: $destRdll not found (some builds place R.dll under bin only)"; }

  Write-Host "Rscript located at: $destRscript"
  & $destRscript --version

  Write-Host "✅ Windows R runtime ready under $destR"
}
finally {
  if (Test-Path $tmp) {
    Write-Host "Cleaning up temp dir $tmp ..."
    Remove-Item -Recurse -Force $tmp
  }
}
