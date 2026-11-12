import http from "node:http";
import fs from "node:fs";
import path from "node:path";
import { spawn } from "node:child_process";

const argv = process.argv.slice(2);

function readFlag(flag, fallback = null) {
  const idx = argv.indexOf(flag);
  if (idx < 0 || idx + 1 >= argv.length) {
    return fallback;
  }
  return argv[idx + 1];
}

const host = readFlag("--host", "127.0.0.1");
const rawPort = readFlag("--port", "3000");
const dataFile = readFlag("--data-file");
const distDir = readFlag("--dist-dir");
const openBrowser = argv.includes("--open");

if (!dataFile || !distDir) {
  console.error("Usage: node tools/report_web_server.mjs --data-file FILE --dist-dir DIR [--host HOST] [--port PORT]");
  process.exit(1);
}

if (!/^(0|[1-9][0-9]{0,4})$/.test(rawPort)) {
  console.error(`[REPORT-WEB] invalid --port value: ${rawPort}`);
  process.exit(1);
}

const port = Number.parseInt(rawPort, 10);
if (port < 0 || port > 65535) {
  console.error(`[REPORT-WEB] invalid --port value: ${rawPort}`);
  process.exit(1);
}

const resolvedDistDir = path.resolve(distDir);
const resolvedDataFile = path.resolve(dataFile);
const resolvedResultsDir = path.dirname(path.dirname(resolvedDataFile));

if (!fs.existsSync(resolvedDataFile)) {
  console.error(`[REPORT-WEB] report data file not found: ${resolvedDataFile}`);
  process.exit(1);
}
if (!fs.existsSync(resolvedDistDir) || !fs.statSync(resolvedDistDir).isDirectory()) {
  console.error(`[REPORT-WEB] web dist directory not found: ${resolvedDistDir}`);
  process.exit(1);
}

const contentType = (filePath) => {
  if (filePath.endsWith(".js")) return "text/javascript; charset=utf-8";
  if (filePath.endsWith(".css")) return "text/css; charset=utf-8";
  if (filePath.endsWith(".json")) return "application/json; charset=utf-8";
  if (filePath.endsWith(".svg")) return "image/svg+xml";
  if (filePath.endsWith(".png")) return "image/png";
  if (filePath.endsWith(".pdf")) return "application/pdf";
  return "text/html; charset=utf-8";
};

const openUrl = (url) => {
  const platform = process.platform;
  const command = platform === "darwin" ? "open" : platform === "win32" ? "cmd" : "xdg-open";
  const args = platform === "win32" ? ["/c", "start", "", url] : [url];
  const child = spawn(command, args, { detached: true, stdio: "ignore" });
  child.on("error", (err) => {
    console.error(`[REPORT-WEB] auto-open failed: ${err.message}`);
    console.error(`[REPORT-WEB] open manually: ${url}`);
  });
  child.unref();
};

const serveFile = (res, filePath) => {
  res.writeHead(200, { "Content-Type": contentType(filePath) });
  fs.createReadStream(filePath).pipe(res);
};

const server = http.createServer((req, res) => {
  const rawPath = (req.url || "/").split("?")[0];
  if (rawPath === "/api/report-data") {
    serveFile(res, resolvedDataFile);
    return;
  }

  if (rawPath.startsWith("/plots/") || rawPath.startsWith("/plot_data/")) {
    const normalizedAssetPath = path.normalize(rawPath).replace(/^(\.\.[/\\])+/, "");
    const candidateAsset = path.join(resolvedResultsDir, normalizedAssetPath);
    if (candidateAsset.startsWith(resolvedResultsDir) && fs.existsSync(candidateAsset) && fs.statSync(candidateAsset).isFile()) {
      serveFile(res, candidateAsset);
      return;
    }
    res.writeHead(404, { "Content-Type": "text/plain; charset=utf-8" });
    res.end("asset not found");
    return;
  }

  const normalizedPath = path.normalize(rawPath).replace(/^(\.\.[/\\])+/, "");
  const requestPath = normalizedPath === "/" ? "/index.html" : normalizedPath;
  const candidate = path.join(resolvedDistDir, requestPath);

  if (candidate.startsWith(resolvedDistDir) && fs.existsSync(candidate) && fs.statSync(candidate).isFile()) {
    serveFile(res, candidate);
    return;
  }

  serveFile(res, path.join(resolvedDistDir, "index.html"));
});

server.listen(port, host, () => {
  const addr = server.address();
  let url = `http://${host}:${port}`;
  if (!addr || typeof addr === "string") {
    console.log(`[REPORT-WEB] url: ${url}`);
  } else {
    url = `http://${host}:${addr.port}`;
    console.log(`[REPORT-WEB] url: ${url}`);
  }
  if (openBrowser) {
    openUrl(url);
  }
});

server.on("error", (err) => {
  console.error(`[REPORT-WEB] server error: ${err.message}`);
  process.exit(1);
});
