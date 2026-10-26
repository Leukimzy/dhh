"""Match author plotting code: Fig. 6c heatmaps use autosomal variants."""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "raw_data" / "fig6c_filtered_sv_hashes.csv"
ALL = ROOT / "data" / "fig6c_all_sites.csv"
AUTO = ROOT / "data" / "fig6c_autosomal_sites.csv"

autosome_ids = set()
with CACHE.open(newline="", encoding="utf-8-sig") as stream:
    for row in csv.DictReader(stream):
        if row["autosome"] == "1":
            autosome_ids.add(row["variant_hash"])

with ALL.open(newline="", encoding="utf-8-sig") as src, \
        AUTO.open("w", newline="", encoding="utf-8") as dst:
    reader = csv.DictReader(src)
    writer = csv.DictWriter(dst, fieldnames=reader.fieldnames)
    writer.writeheader()
    for row in reader:
        if row["variant_hash"] in autosome_ids:
            writer.writerow(row)
print("saved", AUTO)
