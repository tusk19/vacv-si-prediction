# -*- coding: utf-8 -*-
"""
Глубокое покрытие ChEMBL для НЕПОКРЫТЫХ молекул теста:
  SMILES -> similarity 100% (ловит соли/таутомеры) -> активность против VACV
  -> same-doc пары IC50+CC50 -> SI.

Выход (инкрементально): data/chembl_deep_hits.csv  (index, chembl_id, si, kind, doc)
Запуск: python scripts/chembl_deep_coverage.py
"""
import sys
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
    unc = pd.read_csv(DATA / "uncovered_test.csv")
    print(f"Непокрытых: {len(unc)}", flush=True)
    out_path = DATA / "chembl_deep_hits.csv"
    hits = []

    for n, (_, row) in enumerate(unc.iterrows()):
        idx, smi = int(row["index"]), row["SMILES"]
        # 1. 100% similarity -> chembl molecule id
        d = get(f"{BASE}/similarity/{requests.utils.quote(smi, safe='')}/100.json")
        time.sleep(PAUSE)
        mols = (d or {}).get("molecules", [])
        if not mols:
            continue
        mol_ids = [m["molecule_chembl_id"] for m in mols][:5]

        # 2. их активность против VACV
        acts = get(f"{BASE}/activity.json",
                   {"molecule_chembl_id__in": ",".join(mol_ids),
                    "target_chembl_id": VACV, "limit": 100})
        time.sleep(PAUSE)
        alist = (acts or {}).get("activities", [])
        ics = [a for a in alist if a.get("standard_type") in ("IC50", "EC50")
               and a.get("standard_value") and a.get("standard_units") == "nM"]
        if not ics:
            continue
        a0 = min(ics, key=lambda a: float(a["standard_value"]))
        ic50_uM = float(a0["standard_value"]) / 1000.0
        doc = a0.get("document_chembl_id")

        # 3. CC50 из того же документа
        cc50_uM = None
        if doc:
            d2 = get(f"{BASE}/activity.json",
                     {"document_chembl_id": doc,
                      "molecule_chembl_id": a0["molecule_chembl_id"],
                      "standard_type": "CC50", "limit": 100})
            time.sleep(PAUSE)
            ccs = [a for a in (d2 or {}).get("activities", [])
                   if a.get("standard_value") and a.get("standard_units") == "nM"]
            if ccs:
                cc50_uM = float(ccs[0]["standard_value"]) / 1000.0

        if cc50_uM and ic50_uM > 0:
            hits.append({"index": idx, "chembl_id": a0["molecule_chembl_id"],
                         "si": cc50_uM / ic50_uM, "kind": "paired_doc", "doc": doc})
        else:
            hits.append({"index": idx, "chembl_id": a0["molecule_chembl_id"],
                         "si": None, "kind": "ic50_only", "doc": doc})
        if len(unc) and (n + 1) % 20 == 0:
            print(f"  {n+1}/{len(unc)} | попаданий: {sum(1 for h in hits if h['si'])}",
                  flush=True)
            pd.DataFrame(hits).to_csv(out_path, index=False)

    pd.DataFrame(hits).to_csv(out_path, index=False)
    good = sum(1 for h in hits if h["si"])
    print(f"\nГотово: молекул найдено {len(hits)}, с SI: {good} -> {out_path}")


if __name__ == "__main__":
    main()
