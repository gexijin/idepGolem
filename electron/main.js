// main.js
const { app, BrowserWindow, dialog } = require('electron');
const path = require('path');
const fs = require('fs');
const { spawn } = require('child_process');
const net = require('net');
const { fetch } = require('undici');

let childProc = null;

// ---------- logging ----------
const LOG_FILE = path.join(app.getPath('temp'), 'idep-electron.log');
function log(...args) {
  try {
    const line = args.map(x => (typeof x === 'string' ? x : JSON.stringify(x))).join(' ');
    fs.appendFileSync(LOG_FILE, line + '\n');
    console.log(line);
  } catch {}
}

// ---------- crash guards ----------
process.on('uncaughtException', (err) => {
  const msg = (err && err.stack) ? err.stack : String(err);
  log('[uncaughtException]', msg);
  try { dialog.showErrorBox('Uncaught Exception', msg); } catch {}
});
process.on('unhandledRejection', (reason) => {
  const msg = (reason && reason.stack) ? reason.stack : String(reason);
  log('[unhandledRejection]', msg);
});

// ---------- single instance ----------
const gotLock = app.requestSingleInstanceLock();
if (!gotLock) app.quit();
else app.on('second-instance', () => {
  if (global.win) {
    if (global.win.isMinimized()) global.win.restore();
    global.win.focus();
  }
});

// ---------- helpers ----------
function getRuntime() {
  const rp = app.isPackaged ? process.resourcesPath : __dirname;

  if (process.platform === 'win32') {
    const roots = [
      path.join(rp, 'runtime', 'R.win'),
      path.join(rp, 'R.win'),
      path.join(rp, 'resources', 'R.win'),
      path.join(rp, 'runtime', 'R-Portable'),
      path.join(rp, 'R-Portable'),
    ];
    const R_ROOT = roots.find(fs.existsSync);
    const binDir = R_ROOT ? path.join(R_ROOT, 'bin') : null;
    const candidates = [
      binDir && path.join(binDir, 'Rscript.exe'),
      R_ROOT && path.join(R_ROOT, 'bin', 'x64', 'Rscript.exe'),
    ].filter(Boolean);
    const rscript = candidates.find(p => fs.existsSync(p));
    if (!R_ROOT || !rscript) {
      const msg = `Could not locate bundled Rscript.exe.\nresourcesPath: ${rp}\nChecked:\n${roots.join('\n')}\n`;
      log('[FATAL]', msg);
      try { dialog.showErrorBox('Rscript.exe Not Found', msg + `\nLog: ${LOG_FILE}`); } catch {}
      return null;
    }
    const libDir = path.join(R_ROOT, 'library');
    log('[R runtime]', 'R_ROOT=', R_ROOT, 'libDir=', libDir, 'rscript=', rscript);

    return {
      rscript,
      env: {
        R_HOME: R_ROOT,
        R_LIBS: libDir,
        R_LIBS_USER: libDir,
        R_LIBS_SITE: '',
        R_USER: R_ROOT,
        PATH: [binDir, process.env.PATH || ''].filter(Boolean).join(';'),
      },
    };
  }

  if (process.platform === 'darwin') {
    const candidates = [
      path.join(rp, 'runtime', 'R.framework', 'Resources', 'bin', 'Rscript'),
      path.join(rp, 'R.framework', 'Resources', 'bin', 'Rscript'),
      path.join(rp, 'resources', 'R.framework', 'Resources', 'bin', 'Rscript'),
    ];
    const rscript = candidates.find(fs.existsSync);
    if (!rscript) {
      log('[macOS] Rscript not found. Checked:\n' + candidates.join('\n'));
      try { dialog.showErrorBox('Rscript Not Found', 'Bundle R.framework under runtime/.\nSee log: ' + LOG_FILE); } catch {}
      return null;
    }
    const R_RES = path.dirname(path.dirname(rscript)); // .../R.framework/Resources
    return {
      rscript,
      env: {
        R_HOME: R_RES,
        DYLD_FALLBACK_LIBRARY_PATH: path.join(R_RES, 'lib'),
        PATH: [path.join(R_RES, 'bin'), process.env.PATH || ''].filter(Boolean).join(':'),
      },
    };
  }

  // linux
  const roots = [
    path.join(rp, 'runtime', 'R.linux'),
    path.join(rp, 'R.linux'),
    path.join(rp, 'resources', 'R.linux'),
  ];
  const R_ROOT = roots.find(fs.existsSync);
  const binDir = R_ROOT ? path.join(R_ROOT, 'bin') : null;
const rscript = (binDir && fs.existsSync(path.join(binDir, 'Rscript')))
  ? path.join(binDir, 'Rscript')
  : 'Rscript';
  return {
    rscript,
    env: R_ROOT ? {
      R_HOME: R_ROOT,
      LD_LIBRARY_PATH: [path.join(R_ROOT, 'lib'), process.env.LD_LIBRARY_PATH || ''].filter(Boolean).join(':'),
      PATH: [binDir, process.env.PATH || ''].filter(Boolean).join(':'),
    } : {},
  };
}

function isWritableDir(p) {
  try { fs.accessSync(p, fs.constants.W_OK); return fs.statSync(p).isDirectory(); } catch { return false; }
}
function safeKill(proc) {
  try {
    if (proc && !proc.killed) {
      if (process.platform === 'win32') proc.kill('SIGTERM');
      else proc.kill('SIGINT');
    }
  } catch {}
}

async function waitForHttp(url, { timeoutMs = 120000, intervalMs = 500 } = {}) {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    try {
      const ctrl = new AbortController();
      const to = setTimeout(() => ctrl.abort(), Math.min(5000, intervalMs * 4));
      const res = await fetch(url, { method: 'GET', signal: ctrl.signal });
      clearTimeout(to);
      if (res.ok || res.status === 404 || res.status === 403) return true;
    } catch {}
    await new Promise(r => setTimeout(r, intervalMs));
  }
  throw new Error(`Timeout waiting for ${url}`);
}

function getFreePort(start = 7777, end = 7999) {
  return new Promise((resolve, reject) => {
    const server = net.createServer();
    server.unref();
    server.on('error', () => {
      if (start < end) resolve(getFreePort(start + 1, end));
      else reject(new Error('No free ports'));
    });
    server.listen(start, '127.0.0.1', () => {
      const port = server.address().port;
      server.close(() => resolve(port));
    });
  });
}

// update splash progress bar + taskbar progress
function setSplashProgress(progress, statusText) {
  if (!global.win) return;
  try {
    // Taskbar / dock progress
    if (typeof progress === 'number' && progress >= 0 && progress <= 1) {
      global.win.setProgressBar(progress);
    } else {
      global.win.setProgressBar(-1); // clear
    }

    // In-window bar + text
    const pct = typeof progress === 'number' ? Math.round(Math.max(0, Math.min(1, progress)) * 100) : 0;
    const js = `
      if (window.updateSplash) {
        window.updateSplash(${pct}, ${statusText ? JSON.stringify(statusText) : 'null'});
      }
    `;
    global.win.webContents.executeJavaScript(js).catch(() => {});
  } catch {}
}

function showPlaceholder() {
  if (global.win) return;
  global.win = new BrowserWindow({
    width: 900,
    height: 500,
    show: true,
    resizable: true,
    webPreferences: { contextIsolation: true, nodeIntegration: false },
  });

  const html = `
    <html>
    <head>
      <meta charset="utf-8" />
      <title>iDEP – Starting…</title>
    </head>
    <body style="margin:0;padding:0;font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;background:#0f172a;color:#e5e7eb;display:flex;align-items:center;justify-content:center;height:100vh;">
      <div style="max-width:480px;width:100%;padding:24px;">
        <h1 style="margin:0 0 8px 0;font-size:28px;font-weight:600;color:#e5e7eb;">Starting iDEP…</h1>
        <p style="margin:0 0 16px 0;font-size:14px;color:#9ca3af;">
          Initializing embedded R and Shiny server. This may take a moment on first launch.
        </p>

        <div style="margin:16px 0 6px 0;width:100%;height:10px;border-radius:999px;background:#1f2937;overflow:hidden;">
          <div id="bar" style="height:100%;width:8%;background:#22c55e;transition:width 0.35s ease;"></div>
        </div>
        <div style="font-size:12px;color:#9ca3af;">
          <span id="status">Preparing runtime…</span>
        </div>

        <div style="margin-top:16px;font-size:11px;color:#6b7280;line-height:1.4;">
          Log file: <code style="font-size:11px;color:#9ca3af;">${LOG_FILE.replace(/\\/g, '/')}</code><br/>
          If this screen remains for a long time, please open the log file and share it with the developer.
        </div>
      </div>

      <script>
        window.updateSplash = function(pct, text) {
          try {
            var bar = document.getElementById('bar');
            var label = document.getElementById('status');
            if (bar && typeof pct === 'number') {
              var clamped = Math.max(0, Math.min(100, pct));
              bar.style.width = clamped + '%';
            }
            if (label && text) {
              label.textContent = text;
            }
          } catch (e) {
            // ignore
          }
        };
      </script>
    </body>
    </html>
  `;

  global.win.loadURL('data:text/html;charset=utf-8,' + encodeURIComponent(html));
}

// ---------- bootstrap ----------
async function createWindow() {
  const host = '127.0.0.1';
  const port = await getFreePort();
  let shinyPortFromLog = null; // track port reported by Shiny

  const RESOURCES_DIR = process.resourcesPath;
  const APP_DIR = app.isPackaged
    ? path.join(RESOURCES_DIR, 'app')
    : path.join(__dirname, 'app');

  // show splash early
  showPlaceholder();
  setSplashProgress(0.1, 'Preparing data directory…');

  // demo data directory under app
  const DEMO_DIR = path.join(APP_DIR, 'data113');
  let demoDirExists = false;
  try {
    demoDirExists = fs.existsSync(DEMO_DIR) && fs.readdirSync(DEMO_DIR).length > 0;
  } catch {}
  log('[demo data]', 'DEMO_DIR =', DEMO_DIR, 'exists =', demoDirExists);

  // data dir
  const LAUNCH_DIR = process.cwd();
  const overrideDir = process.env.IDEP_DATA_DIR || process.env.IDEP_DATABASE;
  let DATA_PARENT;
  if (overrideDir) DATA_PARENT = path.resolve(overrideDir);
  else if (LAUNCH_DIR && LAUNCH_DIR !== '/' && isWritableDir(LAUNCH_DIR)) DATA_PARENT = path.join(LAUNCH_DIR, 'idep');
  else DATA_PARENT = path.join(app.getPath('userData'), 'idep');
  try { fs.mkdirSync(DATA_PARENT, { recursive: true }); } catch {}

  // sanity
  const appR = path.join(APP_DIR, 'app.R');
  if (!fs.existsSync(appR)) {
    const msg = `Missing app/app.R.\nLooked at: ${appR}\nLog: ${LOG_FILE}`;
    log('[FATAL]', msg);
    try { dialog.showErrorBox('Missing app.R', msg); } catch {}
    app.quit(); return;
  }

  // runtime
  const runtime = getRuntime();
  if (!runtime) { app.quit(); return; }
  const { rscript, env } = runtime;
  setSplashProgress(0.25, 'R runtime located…');

  // OPTIONAL: one-time diagnostics
  try {
    const diag = spawn(rscript, ['-e',
      "cat('LIBPATHS:\\n', paste(.libPaths(), collapse='\\n'), '\\n'); " +
      "cat('ottoPlots available? ', requireNamespace('ottoPlots', quietly=TRUE), '\\n')"
    ], {
      env: { ...process.env, ...env },
      windowsHide: true,
    });
    diag.stdout.on('data', d => log('[R diag stdout]', String(d).trim()));
    diag.stderr.on('data', d => log('[R diag stderr]', String(d).trim()));
  } catch (e) {
    log('[R diag error]', e && e.stack ? e.stack : String(e));
  }

  // write bootstrap
  const bootstrapPath = path.join(DATA_PARENT, 'electron_bootstrap.R');
  const bootstrapSrc = `
args <- commandArgs(trailingOnly = TRUE)
data_dir <- Sys.getenv("IDEP_DATA_DIR", unset = getwd())
app_dir  <- Sys.getenv("IDEP_APP_DIR",  unset = "${APP_DIR.replace(/\\/g, '/')}")
lib_dir  <- Sys.getenv("R_LIBS_USER",   unset = file.path(Sys.getenv("R_HOME","."),"library"))
host     <- Sys.getenv("IDEP_HOST", unset = "127.0.0.1")
port     <- as.integer(Sys.getenv("IDEP_PORT", unset = "${port}"))
demo_dir_hint <- Sys.getenv("IDEP_DEMO_DIR", unset = file.path(app_dir, "data113"))

.libPaths(unique(c(normalizePath(lib_dir, winslash="/", mustWork=FALSE), .libPaths())))
options(shiny.launch.browser = FALSE, golem.app.prod = TRUE)

# ------------------------------------------------------------------------------
# Avoid re-downloading the large demo tarball every launch.
# We cache it under data_dir and reuse it if available.
# ------------------------------------------------------------------------------
cache_tar <- file.path(data_dir, "data113_cache.tar.gz")
demo_url  <- "http://bioinformatics.sdstate.edu/data/data113/data113.tar.gz"

try({
  ns_utils <- asNamespace("utils")
  orig_download_file <- get("download.file", envir = ns_utils)

  patched_download <- function(url, destfile, ...) {
    if (identical(url, demo_url) && file.exists(cache_tar)) {
      message("[bootstrap] Using cached demo tarball: ", cache_tar)
      dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
      file.copy(cache_tar, destfile, overwrite = TRUE)
      return(0L)
    }

    status <- orig_download_file(url, destfile, ...)

    if (identical(url, demo_url) && identical(status, 0L) && file.exists(destfile)) {
      dir.create(dirname(cache_tar), recursive = TRUE, showWarnings = FALSE)
      if (file.copy(destfile, cache_tar, overwrite = TRUE)) {
        message("[bootstrap] Cached demo tarball at: ", cache_tar)
      }
    }

    status
  }

  unlockBinding("download.file", ns_utils)
  assign("download.file", patched_download, envir = ns_utils)
  lockBinding("download.file", ns_utils)
}, silent = TRUE)

# ------------------------------------------------------------------------------
# Avoid re-untarring the demo tarball if data113 already exists and has content.
# ------------------------------------------------------------------------------
try({
  ns_utils2 <- asNamespace("utils")
  orig_untar <- get("untar", envir = ns_utils2)

  patched_untar <- function(tarfile, files = NULL, list = FALSE,
                            exdir = ".", compressed = NA,
                            extras = NULL, verbose = getOption("verbose")) {

    is_demo_tar <- grepl("data113", basename(tarfile), fixed = TRUE)

    if (is_demo_tar &&
        dir.exists(exdir) &&
        length(list.files(exdir, recursive = TRUE)) > 0L) {
      message("[bootstrap] Skipping untar for demo data; existing files in: ", exdir)
      return(invisible(character()))
    }

    orig_untar(tarfile,
               files = files,
               list  = list,
               exdir = exdir,
               compressed = compressed,
               extras = extras,
               verbose = verbose)
  }

  unlockBinding("untar", ns_utils2)
  assign("untar", patched_untar, envir = ns_utils2)
  lockBinding("untar", ns_utils2)
}, silent = TRUE)

# Optional: still log whether we see a pre-extracted demo_dir
if (dir.exists(demo_dir_hint) && length(list.files(demo_dir_hint, recursive = TRUE)) > 0L) {
  options(idep.demo_data_dir = demo_dir_hint)
  message("[bootstrap] Existing demo data dir detected at: ", demo_dir_hint)
} else {
  message("[bootstrap] No existing demo data dir at: ", demo_dir_hint)
}

# log
logfile <- file.path(data_dir, "electron_r.log")
try({
  zz <- file(logfile, open = "a+", encoding = "UTF-8")
  sink(zz, type = "output", split = TRUE)
  sink(zz, type = "message", split = TRUE)
}, silent = TRUE)

# inform Electron about port
writeLines(as.character(port), file.path(data_dir, "idep_port.txt"))

setwd(app_dir)
ok <- TRUE
startup_t0 <- Sys.time()

tryCatch({
  pkg_t0 <- Sys.time()
  if (!requireNamespace("ottoPlots", quietly = TRUE)) {
    stop("Package 'ottoPlots' not found in vendored library: ", paste(.libPaths(), collapse=" | "))
  }
  if (!requireNamespace("idepGolem", quietly = TRUE)) {
    stop("Package 'idepGolem' not found in vendored library: ", paste(.libPaths(), collapse=" | "))
  }
  pkg_t1 <- Sys.time()
  message("[bootstrap] Package load time: ", round(as.numeric(difftime(pkg_t1, pkg_t0, units = "secs")), 2), " s")

  app_t0 <- Sys.time()
  app <- idepGolem::run_app()
  if (!inherits(app, "shiny.appobj")) stop("run_app() did not return a shiny.appobj")
  app_t1 <- Sys.time()
  message("[bootstrap] run_app() time: ", round(as.numeric(difftime(app_t1, app_t0, units = "secs")), 2), " s")

  shiny_t0 <- Sys.time()
  shiny::runApp(app, host = host, port = port, launch.browser = FALSE)
  shiny_t1 <- Sys.time()
  message("[bootstrap] shiny::runApp() returned after ",
          round(as.numeric(difftime(shiny_t1, shiny_t0, units = "secs")), 2), " s")
}, error = function(e) {
  ok <<- FALSE
  cat("FATAL:", conditionMessage(e), "\\n")
}, finally = {
  startup_t1 <- Sys.time()
  message("[bootstrap] Total R-side startup time: ",
          round(as.numeric(difftime(startup_t1, startup_t0, units = "secs")), 2), " s")
  try({ sink(type = "message"); sink(type = "output") }, silent = TRUE)
})

quit(status = if (ok) 0L else 1L, save = "no")
`;
  try { fs.writeFileSync(bootstrapPath, bootstrapSrc, 'utf8'); }
  catch (e) {
    const msg = 'Failed to write bootstrap: ' + (e.stack || String(e));
    log('[bootstrap write error]', msg);
    try { dialog.showErrorBox('R Launch Error', msg + `\n\nLog: ${LOG_FILE}`); } catch {}
    app.quit(); return;
  }

  log(`=== Launch ${new Date().toISOString()} ===`);
  log(`resourcesPath = ${RESOURCES_DIR}`);
  log(`APP_DIR       = ${APP_DIR}`);
  log(`DATA_PARENT   = ${DATA_PARENT}`);
  log(`Rscript       = ${rscript}`);
  log(`bootstrap.R   = ${bootstrapPath}`);

  setSplashProgress(0.35, 'Writing R bootstrap script…');

  // spawn R
  try {
    childProc = spawn(rscript, ['--vanilla', bootstrapPath], {
      cwd: DATA_PARENT,
      env: {
        ...process.env,
        ...env,
        IDEP_DATABASE: DATA_PARENT,
        IDEP_DATA_DIR: DATA_PARENT,
        IDEP_APP_DIR: APP_DIR,
        IDEP_HOST: host,
        IDEP_PORT: String(port),
        IDEP_DEMO_DIR: DEMO_DIR, // pass demo dir hint to R
        R_LIBS_USER: env?.R_LIBS || path.join(path.dirname(rscript), '..', 'library'),
      },
      windowsHide: true,
    });
  } catch (e) {
    const msg = `Failed to spawn Rscript: ${e && e.stack ? e.stack : String(e)}\nLog: ${LOG_FILE}`;
    log('[spawn error]', msg);
    try { dialog.showErrorBox('R Launch Error', msg); } catch {}
    app.quit(); return;
  }

  setSplashProgress(0.5, 'Starting embedded R session…');

  if (childProc && childProc.stdout) {
    childProc.stdout.on('data', d => log('[R stdout]', String(d).trim()));
  }

  if (childProc && childProc.stderr) {
    childProc.stderr.on('data', d => {
      const text = String(d);
      log('[R stderr]', text.trim());

      const m = text.match(/Listening on http:\/\/[^:]+:(\d+)/);
      if (m) {
        shinyPortFromLog = Number(m[1]);
        log(`[port detect] Shiny reports listening on port ${shinyPortFromLog}`);
      }
    });
  }

  childProc.on('close', (code, sig) => {
    log('[R exit]', `code=${code||0}`, sig ? `sig=${sig}` : '');
    if (!app.isQuitting) {
      const html = `
        <html><body style="font-family:sans-serif;padding:16px">
          <h2>Server terminated</h2>
          <p>R exited with code: <b>${code ?? 0}</b> ${sig ? `(signal: ${sig})` : ''}</p>
          <p>See log:</p>
          <pre style="white-space:pre-wrap">${LOG_FILE.replace(/\\/g,'/')}</pre>
        </body></html>`;
      if (!global.win) {
        global.win = new BrowserWindow({
          width: 1200, height: 800, show: true,
          webPreferences: { contextIsolation: true, nodeIntegration: false }
        });
      }
      global.win.loadURL('data:text/html;charset=utf-8,' + encodeURIComponent(html));
    }
  });

  // honor dynamic port file
  const portFile = path.join(DATA_PARENT, 'idep_port.txt');
  let targetPort = port;
  for (let i = 0; i < 60; i++) {
    if (fs.existsSync(portFile)) {
      try {
        const val = fs.readFileSync(portFile, 'utf8').trim();
        if (/^\d+$/.test(val)) {
          targetPort = Number(val);
          log(`Detected dynamic port from R: ${targetPort}`);
        }
        break;
      } catch (e) { log('[portFile read error]', e && e.stack ? e.stack : String(e)); }
    }
    await new Promise(r => setTimeout(r, 500));
  }

  // Prefer the port reported by Shiny
  if (shinyPortFromLog !== null && !Number.isNaN(shinyPortFromLog)) {
    if (shinyPortFromLog !== targetPort) {
    log(`[port override] Using port from Shiny log (${shinyPortFromLog}) instead of port file (${targetPort})`);
    }
    targetPort = shinyPortFromLog;
  }

  const finalURL = `http://${host}:${targetPort}`;
  log(`Final targetURL = ${finalURL}`);
  setSplashProgress(0.7, 'Waiting for Shiny server…');

  try {
    await waitForHttp(finalURL, { timeoutMs: 120000, intervalMs: 500 });
  } catch (err) {
    log('[waitForHttp] Timeout/Error:', err && (err.stack || String(err)));
    try { dialog.showErrorBox('Startup Timeout', `App did not start at ${finalURL} within 120s.\nSee log: ${LOG_FILE}`); } catch {}
    safeKill(childProc);
    return;
  }

  setSplashProgress(0.9, 'Loading user interface…');

  // Show app
  try {
    if (!global.win) {
      global.win = new BrowserWindow({
        width: 1200, height: 800,
        webPreferences: { contextIsolation: true, nodeIntegration: false }
      });
    }
    await global.win.loadURL(finalURL);
    setSplashProgress(-1, ''); // clear taskbar progress
  } catch (e) {
    const msg = `Failed to load ${finalURL}: ${e && e.stack ? e.stack : String(e)}`;
    log('[loadURL error]', msg);
    try { dialog.showErrorBox('Load Error', msg + `\n\nLog: ${LOG_FILE}`); } catch {}
  }

  app.on('before-quit', () => { app.isQuitting = true; safeKill(childProc); });
  app.on('window-all-closed', () => app.quit());
}

app.whenReady().then(createWindow);
