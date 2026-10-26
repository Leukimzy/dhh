"""Join official PanGenie statistics and filters for Nature Fig. 6c."""
from collections import Counter
from pathlib import Path
import csv
import gzip
import hashlib
import math

csv.field_size_limit(100_000_000)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "raw_data"
OUTPUT = ROOT / "data" / "fig6c_all_sites.csv"
FILTERS = SOURCE / "bi_all_filters.tsv.gz"
SUMMARY = SOURCE / "summary_bi_all.tsv.gz"
SELECTED = SOURCE / "fig6c_filtered_sv_hashes.csv"
CLASS_COLUMNS = {
    "SV deletions": "is_large_deletion",
    "SV insertions": "is_large_insertion",
    "SV other": "is_large_complex",
}


def main():
    wanted = {}
    autosomes = set()
    all_counts = Counter()
    auto_counts = Counter()
    if SELECTED.exists():
        with SELECTED.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                variant_hash, category = row["variant_hash"], row["type"]
                wanted[variant_hash] = category
                all_counts[category] += 1
                if row["autosome"] == "1":
                    auto_counts[category] += 1
                    autosomes.add(variant_hash)
    else:
        with gzip.open(FILTERS, "rt", encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream, delimiter="\t")
            head = next(reader)
            ix = {name: head.index(name) for name in ("variant_id", "confidence_level", *CLASS_COLUMNS.values())}
            for i, row in enumerate(reader, 1):
                if int(row[ix["confidence_level"]]) < 1:
                    continue
                category = next((name for name, col in CLASS_COLUMNS.items() if row[ix[col]] == "True"), None)
                if category is None:
                    continue
                variant_id = row[ix["variant_id"]]
                variant_hash = hashlib.blake2b(variant_id.encode(), digest_size=16).hexdigest()
                wanted[variant_hash] = category
                all_counts[category] += 1
                if "chrX" not in variant_id and "chrY" not in variant_id:
                    auto_counts[category] += 1
                    autosomes.add(variant_hash)
                if i % 5_000_000 == 0:
                    print("FILTER ROWS", i, flush=True)
        with SELECTED.open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.writer(stream)
            writer.writerow(("variant_hash", "type", "autosome"))
            writer.writerows((key, category, "1" if key in autosomes else "0")
                             for key, category in wanted.items())
    print("FILTER COUNTS all", dict(all_counts), flush=True)
    print("FILTER COUNTS autosomes", dict(auto_counts), flush=True)

    if not SUMMARY.exists():
        print("SUMMARY NOT YET AVAILABLE; FILTERED IDS SAVED", SELECTED, flush=True)
        return

    matched = Counter()
    missing_af = Counter()
    OUTPUT.parent.mkdir(exist_ok=True)
    with gzip.open(SUMMARY, "rb") as stream, \
            OUTPUT.open("w", encoding="utf-8", newline="") as out:
        head = stream.readline().rstrip(b"\r\n").split(b"\t")
        ix_id = head.index(b"variant_id")
        ix_x = head.index(b"panel_allele_freq")
        ix_y = head.index(b"pangenie-unrelated_allele_freq")
        max_ix = max(ix_id, ix_x, ix_y)
        writer = csv.writer(out)
        writer.writerow(("variant_hash", "type", "AF_assemblies", "AF_PanGenie"))
        for i, line in enumerate(stream, 1):
            if i % 5_000_000 == 0:
                print("SUMMARY ROWS", i, flush=True)
            row = line.split(b"\t", max_ix+1)
            variant_id = row[ix_id]
            if b"-SNV-" in variant_id:
                continue
            variant_hash = hashlib.blake2b(variant_id, digest_size=16).hexdigest()
            category = wanted.get(variant_hash)
            if category is None:
                continue
            try:
                x = float(row[ix_x]); y = float(row[ix_y])
            except (ValueError, IndexError):
                missing_af[category] += 1
                continue
            if not (math.isfinite(x) and math.isfinite(y)):
                missing_af[category] += 1
                continue
            writer.writerow((variant_hash, category, x, y))
            matched[category] += 1
    print("MATCHED", dict(matched), flush=True)
    print("MISSING AF", dict(missing_af), flush=True)
    print("SAVED", OUTPUT, flush=True)


if __name__ == "__main__":
    main()
