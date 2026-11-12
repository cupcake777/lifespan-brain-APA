# Export the three summary sheets used by ED6 and ED9 from the supplementary workbook to TSV.
import csv, os
import openpyxl

ROOT = os.environ["APA_ROOT"]
xlsx = f"{ROOT}/data/04_qtl/Supplementary_Tables_draft.xlsx"
wanted = {"S_replication_pi1": "ed6", "S_cell_composition": "ed6", "S_MAGMA": "ed9"}
wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
for sheet, fig in wanted.items():
    ws = wb[sheet]
    path = os.path.join(ROOT, "data/figures", fig, f"{sheet}.tsv")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for row in ws.iter_rows(values_only=True):
            w.writerow(["" if v is None else v for v in row])
