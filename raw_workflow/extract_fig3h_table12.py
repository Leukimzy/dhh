"""Extract Fig. 3h values from the paper's Supplementary Table 12.

The main figure has five legend classes. The source table additionally lists a
small 'Low confidence' class (0.03--0.05%); group it with 'Missing' for display
and retain its original counts in a separate six-class audit table.
"""
from pathlib import Path
import csv
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "raw_data" / "nature_supp_tables.xlsx"
OUT = ROOT / "data"
METHOD_COLUMNS = {
    "HiFi-DeepVariant": (1, 2),
    "Dipcall": (3, 4),
    "MC": (5, 6),
    "PGGB": (7, 8),
}
SOURCE_CLASSES = {
    2: "Missing",
    3: "Low confidence",
    4: "Singleton",
    5: "Rare",
    6: "Low frequency",
    7: "Common",
}
DISPLAY_CLASSES = ("Common", "Low frequency", "Rare", "Singleton", "Missing")


def main():
    wb = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    ws = wb["Supplementary Table 12"]
    rows = list(ws.iter_rows(min_row=3, max_row=9, min_col=1, max_col=9, values_only=True))
    raw, shown = [], []
    for method, (count_col, pct_col) in METHOD_COLUMNS.items():
        total = int(str(rows[6][count_col]).replace(",", ""))
        by_class = {}
        for row_index, category in SOURCE_CLASSES.items():
            row = rows[row_index - 2]
            count = int(str(row[count_col]).replace(",", ""))
            by_class[category] = count
            raw.append((method, category, count, str(row[pct_col]), total))
        assert sum(by_class.values()) == total, (method, by_class, total)
        by_class["Missing"] += by_class["Low confidence"]
        for category in DISPLAY_CLASSES:
            count = by_class[category]
            shown.append((method, category, count, total, count / total * 100))
    wb.close()

    with (OUT / "fig3h_six_class_counts.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(("method", "category", "count", "source_percent", "total")); w.writerows(raw)
    with (OUT / "fig3h_five_classes.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(("method", "category", "count", "total", "percent")); w.writerows(shown)
    for method in METHOD_COLUMNS:
        print(method, [(r[1], round(r[4], 4)) for r in shown if r[0] == method])


if __name__ == "__main__":
    main()
