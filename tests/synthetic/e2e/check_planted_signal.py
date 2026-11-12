#!/usr/bin/env python3
# Check that the module 02 chain recovers the signal planted by make_synthetic_cohort.py.
# Usage: check_planted_signal.py RUN_DIR truth.json   (reads RUN_DIR/outputs/02_sites.tsv)

import json
import sys

import pandas as pd

# --- inputs ---
run, truth_path = sys.argv[1], sys.argv[2]
truth = json.load(open(truth_path))
sites = pd.read_csv(f"{run}/outputs/02_sites.tsv", sep="\t")
sites["dyn"] = sites["dynamic_pass"].astype(str).str.upper().eq("TRUE")
lengthen = set(truth["lengthening_sites"])
shorten = set(truth["shortening_sites"])
planted = lengthen | shorten
null = sites[~sites.site_id.isin(planted)]
up = sites[sites.site_id.isin(lengthen)]
dn = sites[sites.site_id.isin(shorten)]

# --- rates and assertions ---
# pass needs >= 60% recovery of each planted direction and a null rate below the planted rate; exit 1 otherwise
rate_up = (up.dyn & (up.delta_logit > 0)).mean()
rate_dn = (dn.dyn & (dn.delta_logit < 0)).mean()
rate_null = null.dyn.mean()
print(f"lengthening recovered (dynamic, delta>0): {rate_up:.3f} of {len(up)}")
print(f"shortening recovered (dynamic, delta<0): {rate_dn:.3f} of {len(dn)}")
print(f"dynamic_pass rate among null sites: {rate_null:.3f} of {len(null)}")
ok = rate_up >= 0.6 and rate_dn >= 0.6 and rate_null < min(rate_up, rate_dn)
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
