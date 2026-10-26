"""Count small-variant sites directly from the released Year1 graph VCFs.

Each VCF record is one variant site. This uses the *processed* v1.0
vcfbub/wave callsets, so the counts need not equal the paper's pre-wave
Fig. 3b pipeline exactly. No bar height is estimated from the figure.
"""
from collections import Counter
from pathlib import Path
import csv
import gzip
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "raw_data"
OUTPUT = ROOT / "data" / "fig3b_vcf_counts.csv"
DATA = {
    "MC": SOURCE / "mc_v1_processed.vcf.gz",
    "PGGB": SOURCE / "pggb_v1_processed.vcf.gz",
}


def classify(ref: bytes, alt: bytes):
    alleles = alt.split(b",")
    if any(a in (b".", b"*") or a.startswith(b"<") for a in alleles):
        return None
    lengths = [len(ref), *(len(a) for a in alleles)]
    if max(lengths) >= 50:
        return None
    if all(n == 1 for n in lengths):
        return "SNP"
    if len(set(lengths)) == 1:
        return "MNP"
    return "INDEL"


def count_one(path: Path):
    counts = Counter()
    diagnostics = Counter()
    with gzip.open(path, "rb") as stream:
        for i, line in enumerate(stream, 1):
            if line.startswith(b"#"):
                continue
            columns = line.split(b"\t", 8)
            if len(columns) < 5:
                diagnostics["malformed"] += 1
                continue
            diagnostics["records"] += 1
            category = classify(columns[3], columns[4])
            if category is None:
                diagnostics["not_small_or_symbolic"] += 1
                continue
            n_alleles = columns[4].count(b",") + 2
            allele_bin = "2" if n_alleles == 2 else "3" if n_alleles == 3 else ">3"
            counts[category, allele_bin] += 1
            if diagnostics["records"] % 5_000_000 == 0:
                print(path.name, diagnostics["records"], flush=True)
    return counts, diagnostics


def main():
    methods = sys.argv[1:] or list(DATA)
    existing = {}
    if OUTPUT.exists():
        with OUTPUT.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                existing[row["method"], row["type"], row["alleles"]] = int(row["sites"])
    for method in methods:
        if method not in DATA:
            raise ValueError(method)
        counts, diagnostics = count_one(DATA[method])
        print(method, json.dumps({"counts": {f"{k[0]}_{k[1]}":v for k,v in counts.items()}, "diagnostics":diagnostics}), flush=True)
        for category in ("SNP", "MNP", "INDEL"):
            for allele_bin in ("2", "3", ">3"):
                existing[method, category, allele_bin] = counts[category, allele_bin]
        OUTPUT.parent.mkdir(exist_ok=True)
        with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("method", "type", "alleles", "sites"))
            for m in ("MC", "PGGB"):
                for cat in ("SNP", "MNP", "INDEL"):
                    for ab in ("2", "3", ">3"):
                        if (m, cat, ab) in existing:
                            writer.writerow((m, cat, ab, existing[m, cat, ab]))
        print("saved", OUTPUT, flush=True)


if __name__ == "__main__":
    main()
