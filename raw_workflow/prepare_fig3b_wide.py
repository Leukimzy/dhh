"""Turn the 18 raw VCF counts into the six-row table read by reproduce.R."""
from pathlib import Path
import csv

root = Path(__file__).resolve().parents[1]
source = root / "data" / "fig3b_vcf_counts.csv"
target = root / "data" / "fig3b_counts_wide.csv"

with source.open(newline="", encoding="utf-8-sig") as stream:
    values = {
        (r["method"], r["type"], r["alleles"]): int(r["sites"])
        for r in csv.DictReader(stream)
    }

with target.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow(("method", "type", "two", "three", "more"))
    for category in ("SNP", "MNP", "INDEL"):
        for method in ("MC", "PGGB"):
            writer.writerow((method, category,
                             values[method, category, "2"],
                             values[method, category, "3"],
                             values[method, category, ">3"]))
print(target)
