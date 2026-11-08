# -*- coding: utf-8 -*-
"""
Прямой свип Selectivity Index из ChEMBL:
  1. Все молекулы теста -> ChEMBL ID (similarity 100).
  2. Для всех найденных ID -> activity standard_type='Selectivity Index'.
  3. Фильтр доверия: assay_description содержит vaccinia/orthopox/vero/cowpox
     ИЛИ документ из списка VACV-документов.
Выход: data/chembl_si_direct_hits.csv (index, chembl_id, si, assay, doc)
"""
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
BASE = "https://www.ebi.ac.uk/chembl/api/data"
PAUSE = 0.34
TRUST = ("vaccinia", "orthopox", "vero", "cowpox", "variola", "poxvirus")

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
    te = pd.read_csv(DATA / "test.csv")
    print("== Шаг 1: test -> ChEMBL ID ==", flush=True)
    id_map = {}  # index -> chembl_id
    for n, (_, r) in enumerate(te.iterrows()):
        d = get(f"{BASE}/similarity/{requests.utils.quote(r['SMILES'], safe='')}/100.json")
        time.sleep(PAUSE)
        mols = (d or {}).get("molecules", [])
        if mols:
            id_map[int(r["index"])] = mols[0]["molecule_chembl_id"]
        if (n + 1) % 50 == 0:
            print(f"  {n+1}/{len(te)} | найдено {len(id_map)}", flush=True)
    print(f"ChEMBL ID для {len(id_map)} молекул теста")

    print("== Шаг 2: прямые Selectivity Index ==", flush=True)
    ids = sorted(set(id_map.values()))
    hits = []
    for i in range(0, len(ids), 40):
        chunk = ids[i:i + 40]
        d = get(f"{BASE}/activity.json",
                {"molecule_chembl_id__in": ",".join(chunk),
                 "standard_type": "Selectivity Index", "limit": 1000})
        time.sleep(PAUSE)
        for a in (d or {}).get("activities", []):
            if not a.get("standard_value"):
                continue
            desc = (a.get("assay_description") or "").lower()
            if any(k in desc for k in TRUST):
                hits.append({"chembl_id": a["molecule_chembl_id"],
                             "si": float(a["standard_value"]),
                             "assay": (a.get("assay_description") or "")[:120],
                             "doc": a.get("document_chembl_id")})
        print(f"  {min(i+40, len(ids))}/{len(ids)} | SI-записей: {len(hits)}", flush=True)

    df = pd.DataFrame(hits)
    if not len(df):
        print("Прямых SI не найдено.")
        return
    rev = {v: k for k, v in id_map.items()}
    df["index"] = df["chembl_id"].map(rev)
    df = df.dropna(subset=["index"])
    df["index"] = df["index"].astype(int)
    df = df.groupby(["index", "chembl_id"], as_index=False).agg(
        si=("si", "median"), assay=("assay", "first"), doc=("doc", "first"))
    df.to_csv(DATA / "chembl_si_direct_hits.csv", index=False)
    print(f"\nГотово: {len(df)} молекул теста с прямым SI -> data/chembl_si_direct_hits.csv")
    print(df["si"].describe().to_string())


if __name__ == "__main__":
    main()
