#!/usr/bin/env python3
"""Static checks for the paper companion repository."""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
errors: list[str] = []
warnings: list[str] = []

required = [
    "README.md",
    "docs/figure_map.md",
    "docs/input_layout.md",
    "config/paths.example.sh",
    "scripts/render_figures.sh",
    "figure_generation/shared/theme.R",
    *[f"figure_generation/main_figures/figure_{i}.R" for i in range(1, 7)],
    *[f"figure_generation/extended_data_figures/extended_data_{i}.R" for i in range(1, 10)],
    *[f"analysis/{name}/README.md" for name in (
        "01_apa_quantification",
        "02_lifespan_trajectories",
        "03_multiomics_integration",
        "04_stage_specific_qtl",
        "05_genotype_by_stage",
        "06_disease_genetics",
    )],
]
for rel in required:
    if not (ROOT / rel).is_file():
        errors.append(f"missing required file: {rel}")

# Python syntax.
for path in ROOT.rglob("*.py"):
    try:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except SyntaxError as exc:
        errors.append(f"python syntax: {path.relative_to(ROOT)}:{exc.lineno}: {exc.msg}")

# The rewrite intentionally removes internal snapshot paths and development versions.
for path in ROOT.rglob("*"):
    if not path.is_file() or ".git" in path.parts or path.suffix.lower() in {".png", ".pdf"}:
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        continue
    rel = path.relative_to(ROOT)
    if rel == Path("scripts/check_repository.py"):
        continue
    if re.search(r"/(?:root|gpfs|mnt)/", text):
        errors.append(f"absolute internal path in {rel}")
    if rel.parts[:2] in {
        ("figure_generation", "main_figures"),
        ("figure_generation", "extended_data_figures"),
    }:
        if "release/figures" in text or "APA_WORK_ROOT" in text:
            errors.append(f"legacy figure-tree dependency in {rel}")
    if re.search(r"(?:composite_)?v\d+", path.name, re.I):
        errors.append(f"versioned development filename retained: {rel}")

# No bundled third-party products or generated outputs.
for banned in ("easycoloc", "ldscore", "node_modules"):
    hits = [p for p in ROOT.rglob(banned)]
    if hits:
        errors.append(f"vendored/generated directory retained: {hits[0].relative_to(ROOT)}")

# Every Markdown link to a repository-local file or directory must resolve.
for path in ROOT.rglob("*.md"):
    text = path.read_text(encoding="utf-8")
    for target in re.findall(r"\[[^]]*\]\(([^)]+)\)", text):
        target = target.split("#", 1)[0]
        if not target or re.match(r"(?:https?|mailto):", target):
            continue
        if not (path.parent / target).resolve().exists():
            errors.append(f"broken Markdown link in {path.relative_to(ROOT)}: {target}")

print(f"repository: {ROOT}")
print(f"files: {sum(p.is_file() and '.git' not in p.parts for p in ROOT.rglob('*'))}")
print(f"python syntax errors: {sum(e.startswith('python syntax') for e in errors)}")
for item in warnings:
    print(f"WARNING: {item}")
for item in errors:
    print(f"ERROR: {item}")
if errors:
    sys.exit(1)
print("PASS: static repository checks")
