#!/usr/bin/env python3
# Download POSTAR3 CLIP-seq binding sites for each RBP in the browse list, one table per RBP and CLIP technology.
# Writes per-RBP tables and one merged table to data/07_coloc, and a download log to logs/.

# --- packages ---
import json
import os
import random
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests


# --- config ---
ROOT = Path(os.environ["APA_ROOT"])
BASE_URL = "http://111.198.139.65/script/ajax/RBS_show_CLIP_result.php"  # POSTAR3 mirror
# main site: http://postar.ncrnalab.org/script/ajax/RBS_show_CLIP_result.php
HEADERS = {
    "Accept": "application/json, text/javascript, */*; q=0.01",
    "Accept-Language": "zh-CN,zh;q=0.9,en;q=0.8",
    "Referer": "http://111.198.139.65/RBS.html",
    "Origin": "http://111.198.139.65",
    "X-Requested-With": "XMLHttpRequest",
}
CLIP_TABLES = [
    "human_clipdb_eclip",      # ENCODE eCLIP (preferred, highest quality)
    "human_clipdb_ce1",        # eCLIP PureCLIP
    "human_clipdb_ce2",        # eCLIP Piranha
    "human_clipdb_i1",         # iCLIP CITS
    "human_clipdb_i2",         # iCLIP Piranha
    "human_clipdb_i3",         # iCLIP PureCLIP
    "human_clipdb_par1",       # PAR-CLIP PARalyzer
    "human_clipdb_par2",       # PAR-CLIP Piranha
    "human_clipdb_par3",       # PAR-CLIP MiClip
    "human_clipdb_hits1",      # HITS-CLIP CIMS
    "human_clipdb_hits2",      # HITS-CLIP Piranha
    "human_clipdb_hits3",      # HITS-CLIP CLIPper
    "human_clipdb_hits4",      # HITS-CLIP CTK
]
RBP_CSV = ROOT / "ref/POSTAR3_CLIPdb_module_browse_RBP_sub_info.csv"
SAVE_DIR = ROOT / "ref/postar3_all_RBPs_data"
LOG_FILE = ROOT / "logs/download_log.txt"
RETRY = 3


# --- helpers ---
def polite_delay(base_seconds, variance):
    delay = base_seconds + random.uniform(0, variance)
    time.sleep(delay)


# --- load RBP list ---
rbp_list = pd.read_csv(RBP_CSV).iloc[:, 0].dropna().tolist()
rbp_list = [str(rbp).strip() for rbp in rbp_list]

SAVE_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

all_data = []
success_count = 0
total_records = 0

with LOG_FILE.open("w", encoding="utf-8") as log:
    log.write(f"download start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log.write(f"RBP count: {len(rbp_list)}\n")

    for idx, rbp in enumerate(rbp_list, 1):
        log.write(f"\n[{idx}/{len(rbp_list)}] RBP: {rbp}\n")

        rbp_dir = SAVE_DIR / rbp
        rbp_dir.mkdir(parents=True, exist_ok=True)

        rbp_records = 0
        rbp_success_tables = 0

        for table in CLIP_TABLES:
            tech = table.replace("human_clipdb_", "").upper()
            payload = {"species": "human", "geneName": rbp, "table": table}

            # one request per RBP and CLIP table, with retries
            df = None
            count = 0
            for attempt in range(RETRY):
                try:
                    resp = requests.post(BASE_URL, data=payload, headers=HEADERS, timeout=60)
                    resp.raise_for_status()
                    data = resp.json()
                    if not data:
                        break
                    df = pd.DataFrame(data)
                    df["RBP_Query"] = rbp
                    df["CLIP_Technology"] = tech
                    df["Download_Time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    count = len(df)
                    break
                except requests.exceptions.Timeout:
                    if attempt < RETRY - 1:
                        time.sleep(10 + random.uniform(5, 10))
                except requests.exceptions.RequestException:
                    if attempt < RETRY - 1:
                        time.sleep(5 + random.uniform(3, 7))
                except json.JSONDecodeError:
                    break
                except Exception:
                    break

            if df is not None and count > 0:
                filename = rbp_dir / f"{rbp}_{tech}.csv"
                df.to_csv(filename, index=False, encoding="utf-8-sig")
                log.write(f"    {tech}: {count} records\n")
                all_data.append(df)
                rbp_records += count
                rbp_success_tables += 1
                success_count += 1
            else:
                log.write(f"    {tech}: no data\n")

            # pause between CLIP tables of the same RBP
            polite_delay(5, 3)

        log.write(f"  subtotal: {rbp_records} records ({rbp_success_tables}/{len(CLIP_TABLES)} CLIP tables)\n")
        total_records += rbp_records

        # longer pause between RBPs
        polite_delay(10, 5)

        if idx % 50 == 0:
            log.write(f"\n[progress] {idx}/{len(rbp_list)} - total records: {total_records}\n")

    # merged table
    if all_data:
        merged_df = pd.concat(all_data, ignore_index=True)
        merged_file = SAVE_DIR / "ALL_RBP_binding_sites_merged.csv"
        merged_df.to_csv(merged_file, index=False, encoding="utf-8-sig")
        log.write(f"\nmerged file: {merged_file}\n")
        log.write(f"merged records: {len(merged_df)}\n")
    else:
        log.write("\nWARNING: no data retrieved\n")

    log.write(f"\ndownload finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log.write(f"total RBPs: {len(rbp_list)}\n")
    log.write(f"successful CLIP tables: {success_count}\n")
    log.write(f"total records: {total_records}\n")
