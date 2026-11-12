#!/usr/bin/env python3
# Baseline preflight: required paths and schemas must pass before formal analysis.
# Optional extensions (eQTL coloc, protein metadata, PFC RBP, GMT sets) may be absent.

from __future__ import annotations

import gzip
import importlib.util
import json
import os
from pathlib import Path
from zipfile import ZipFile

import pyarrow.parquet as pq
import yaml

# --- paths ---
ROOT = Path(os.environ["APA_ROOT"])
CONFIG = Path(__file__).resolve().parent / "config.yml"
OUT_DIR = ROOT / "data" / "07_coloc"
LIFESTAGE = ROOT / "data" / "05_gxs" / "LifeStage.txt"
OPTIONAL_EXPECTATIONS = {
    "inputs.eqtl_coloc_results": "pending_extension",
    "inputs.protein_sample_metadata": "not_available",
    "reference_resources.pfc_rbp_expression": "not_available",
    "reference_resources.enrichment_gene_sets": "use_package_resources",
}


# --- helpers ---
def _expand_env(obj):
    """Expand ${VAR} in config values (paths set via config/paths.sh)."""
    if isinstance(obj, dict):
        return {k: _expand_env(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env(v) for v in obj]
    return os.path.expandvars(obj) if isinstance(obj, str) else obj


def path_exists(path_text: str) -> bool:
    if path_text.startswith("python_module:"):
        module = path_text.split(":", 1)[1]
        return importlib.util.find_spec(module) is not None
    return Path(path_text).exists()


def read_header(path_text: str) -> list[str]:
    path = Path(path_text)
    if path_text.endswith(".parquet"):
        return pq.ParquetFile(path_text).schema_arrow.names
    if path_text.endswith(".gz"):
        with gzip.open(path, "rt", errors="replace") as handle:
            line = handle.readline().rstrip("\n\r")
    elif path_text.endswith(".zip"):
        with ZipFile(path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if not members:
                return []
            line = archive.open(members[0]).readline().decode("utf-8", errors="replace").rstrip("\n\r")
    else:
        with path.open("rt", errors="replace") as handle:
            line = handle.readline().rstrip("\n\r")
    delimiter = "\t" if "\t" in line else "," if "," in line else None
    return line.split(delimiter) if delimiter else line.split()


def add(rows: list[dict], section: str, key: str, path: str, status: str, requirement: str) -> None:
    rows.append(
        {
            "Section": section,
            "Key": key,
            "Path": path,
            "Requirement": requirement,
            "Status": status,
        }
    )


# --- config and path checks ---
config = _expand_env(yaml.safe_load(CONFIG.read_text()))
rows: list[dict] = []
failures: list[dict] = []

for section in ("inputs", "reference_resources"):
    for key, meta in (config.get(section) or {}).items():
        if not (isinstance(meta, dict) and "path" in meta):
            continue
        path = str(meta.get("path") or "").strip()
        fq_key = f"{section}.{key}"
        required = bool(meta.get("required_for_baseline", True))
        if not path:
            if required:
                add(rows, section, key, path, "FAIL missing required baseline path", "required")
                failures.append(rows[-1])
            else:
                status = f"OPTIONAL {meta.get('absent_status', OPTIONAL_EXPECTATIONS.get(fq_key, 'not_available'))}"
                add(rows, section, key, path, status, "optional_extension")
            continue
        status = "OK" if path_exists(path) else "FAIL missing path"
        add(rows, section, key, path, status, "required" if required else "optional_extension")
        if required and status.startswith("FAIL"):
            failures.append(rows[-1])

ld_path = str(((config.get("reference_resources") or {}).get("ld_reference") or {}).get("path") or "")
if ld_path.endswith(".pgen"):
    pgen = Path(ld_path)
    missing = [str(pgen.with_suffix(ext)) for ext in (".pvar", ".psam") if not pgen.with_suffix(ext).exists()]
    status = "OK .pgen/.pvar/.psam present" if not missing else "FAIL missing " + ",".join(missing)
    add(rows, "reference_resources", "ld_reference_companions", str(pgen.with_suffix("")), status, "required")
    if status.startswith("FAIL"):
        failures.append(rows[-1])

if LIFESTAGE.exists():
    header = read_header(str(LIFESTAGE))
    status = "OK columns_present sample,LifeStage" if {"sample", "LifeStage"}.issubset(header) else "FAIL missing sample/LifeStage"
else:
    status = "FAIL missing path"
add(rows, "extra_key_input", "LifeStage", str(LIFESTAGE), status, "required")
if status.startswith("FAIL"):
    failures.append(rows[-1])

# --- schema checks ---
for key, meta in (config.get("inputs") or {}).items():
    if not isinstance(meta, dict):
        continue
    required_columns = meta.get("required_columns") or []
    path = str(meta.get("path") or "").strip()
    required = bool(meta.get("required_for_baseline", True))
    if not required_columns:
        continue
    if not path:
        status = "SKIP optional absent path" if not required else "FAIL missing required baseline path"
        add(rows, "schema", key, path, status, "required" if required else "optional_extension")
        if required:
            failures.append(rows[-1])
        continue
    if not path_exists(path):
        status = "SKIP missing optional path" if not required else "FAIL missing required baseline path"
        add(rows, "schema", key, path, status, "required" if required else "optional_extension")
        if required:
            failures.append(rows[-1])
        continue
    header = read_header(path)
    missing = [col for col in required_columns if col not in header]
    status = "OK required_columns_present" if not missing else "FAIL missing_columns " + ",".join(missing)
    status += " | header=" + ",".join(header[:16])
    add(rows, "schema", key, path, status, "required" if required else "optional_extension")
    if required and status.startswith("FAIL"):
        failures.append(rows[-1])

# --- outputs ---
cols = ("Section", "Key", "Path", "Requirement", "Status")
lines = ["\t".join(cols)] + ["\t".join(str(row[col]) for col in cols) for row in rows]
OUT_DIR.mkdir(parents=True, exist_ok=True)
(OUT_DIR / "preflight_status.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
print(f"\nSUMMARY\tfailures\t{len(failures)}")
summary = {
    "failures": len(failures),
    "optional_extensions": {
        row["Key"]: row["Status"]
        for row in rows
        if row["Requirement"] == "optional_extension" and row["Path"] == ""
    },
}
(OUT_DIR / "preflight_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
if failures:
    raise SystemExit(1)
