# -*- coding: utf-8 -*-
"""
Добивка: для молекул из chembl_deep_hits.csv со статусом ic50_only
ищем CC50 по ВСЕМУ ChEMBL (любой документ), SI = CC50/IC50.
IC50 берём из активности против VACV (из первого скрипта — пересчитываем).

Выход: data/chembl_deep_hits2.csv (обновлённый: index, chembl_id, si, kind, doc)
Запуск: python scripts/chembl_deep_cc50.py
"""
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
BASE = "https://www.ebi.ac.uk/chembl/api/data"
VACV = "CHEMBL613493"
PAUSE = 0.34

S = requests.Session()
S.headers.update({"User-Agent": "nuclear-it-hack/1.0 (educational)"})


def get(url, params=None):
    for attempt in range(4):
        try:
            r = S.get(url, params=params, timeout=60)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
        time.sleep(1.5 * (attempt + 1))
    return None


def main():
    hits = pd.read_csv(DATA / "chembl_deep_hits.csv")
    only = hits[hits["kind"] == "ic50_only"]
    print(f"ic50_only молекул: {len(only)}")
    if not len(only):
        return

    for i, row in only.iterrows():
        mid = row["chembl_id"]
        # IC50 против VACV (минимальный)
        a = get(f"{BASE}/activity.json",
                {"molecule_chembl_id": mid, "target_chembl_id": VACV, "limit": 100})
        time.sleep(PAUSE)
        ics = [x for x in (a or {}).get("activities", [])
               if x.get("standard_type") in ("IC50", "EC50")
               and x.get("standard_value") and x.get("standard_units") == "nM"]
        if not ics:
            continue
        ic50 = min(float(x["standard_value"]) for x in ics)
        # CC50 по всему ChEMBL (только точные "=")
        c = get(f"{BASE}/activity.json",
                {"molecule_chembl_id": mid, "standard_type": "CC50", "limit": 100})
        time.sleep(PAUSE)
        ccs = [x for x in (c or {}).get("activities", [])
               if x.get("standard_value") and x.get("standard_units") == "nM"
               and x.get("standard_relation") == "="]
        if not ccs:
            continue
        import statistics
        cc50 = statistics.median(float(x["standard_value"]) for x in ccs)
        if ic50 > 0:
            hits.loc[i, "si"] = cc50 / ic50
            hits.loc[i, "kind"] = "global_cc50"
            print(f"  {mid}: SI = {cc50/ic50:.1f}", flush=True)

    hits.to_csv(DATA / "chembl_deep_hits.csv", index=False)
    good = hits["si"].notna().sum()
    print(f"\nИтого с SI: {good} (обновлён chembl_deep_hits.csv)")


if __name__ == "__main__":
    main()
