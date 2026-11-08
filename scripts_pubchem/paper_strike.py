# -*- coding: utf-8 -*-
"""
Хирургический удар: находим статью train-монстра (SI 36067) в ChEMBL,
выкачиваем ВСЕ измерения из её документов и маппим на тест.

Выход: data/paper_strike_hits.csv (index, chembl_id, si, kind, doc, doi)
"""
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from features import canonicalize  # noqa: E402

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
    tr = pd.read_csv(DATA / "train_v2.csv").dropna(subset=["SI VACV"])
    mon = tr[tr["SI VACV"] > 2000]

    # 1. документы монстров
    docs = set()
    for _, r in mon.iterrows():
        d = get(f"{BASE}/similarity/{requests.utils.quote(r['SMILES'], safe='')}/100.json")
        time.sleep(PAUSE)
        mols = (d or {}).get("molecules", [])
        if not mols:
            print(f"  {r['Title']}: нет в ChEMBL")
            continue
        mid = mols[0]["molecule_chembl_id"]
        acts = get(f"{BASE}/activity.json",
                   {"molecule_chembl_id": mid, "limit": 200})
        time.sleep(PAUSE)
        dcs = {a["document_chembl_id"] for a in (acts or {}).get("activities", [])
               if a.get("document_chembl_id")}
        print(f"  {r['Title']} (SI {r['SI VACV']:.0f}) -> {mid}, документов: {len(dcs)}")
        docs |= dcs
    print(f"\nВсего документов: {len(docs)}")

    # 2. DOI документов
    doi = {}
    for doc in docs:
        d = get(f"{BASE}/document/{doc}.json")
        time.sleep(PAUSE)
        if d:
            doi[doc] = d.get("doi") or ""
    print("DOI:", doi)

    # 3. все активности из этих документов, связанные с VACV/токсичностью
    rows = []
    for doc in docs:
        d = get(f"{BASE}/activity.json", {"document_chembl_id": doc, "limit": 1000})
        time.sleep(PAUSE)
        for a in (d or {}).get("activities", []):
            st = a.get("standard_type")
            if st in ("IC50", "EC50", "CC50", "Selectivity Index") and a.get("standard_value"):
                rows.append({"chembl_id": a["molecule_chembl_id"],
                             "smiles": a.get("canonical_smiles") or "",
                             "type": st, "value": float(a["standard_value"]),
                             "units": a.get("standard_units"), "doc": doc,
                             "doi": doi.get(doc, "")})
    df = pd.DataFrame(rows)
    print(f"\nАктивностей выкачано: {len(df)}")
    if not len(df):
        return

    # 4. SI на молекулу: прямой или пара
    out = []
    for mid, g in df.groupby("chembl_id"):
        si_direct = g[g["type"] == "Selectivity Index"]["value"]
        smi = g["smiles"].iloc[0]
        doc = g["doc"].iloc[0]
        if len(si_direct):
            out.append({"chembl_id": mid, "smiles": smi, "si": float(si_direct.median()),
                        "kind": "direct", "doc": doc, "doi": g["doi"].iloc[0]})
            continue
        ics = g[(g["type"].isin(["IC50", "EC50"])) & (g["units"] == "nM")]["value"]
        ccs = g[(g["type"] == "CC50") & (g["units"] == "nM")]["value"]
        if len(ics) and len(ccs) and ics.median() > 0:
            out.append({"chembl_id": mid, "smiles": smi,
                        "si": float(ccs.median() / ics.median()), "kind": "paired",
                        "doc": doc, "doi": g["doi"].iloc[0]})
    si_df = pd.DataFrame(out)
    print(f"Молекул с SI в статьях монстров: {len(si_df)}")
    print(si_df["si"].describe().to_string())

    # 5. маппинг на тест (канонический SMILES)
    te = pd.read_csv(DATA / "test.csv")
    te["canon"] = te["SMILES"].map(canonicalize)
    si_df["canon"] = si_df["smiles"].map(canonicalize)
    m = te.merge(si_df, on="canon", how="inner")
    m = m[["index", "Title", "chembl_id", "si", "kind", "doc", "doi"]]
    m.to_csv(DATA / "paper_strike_hits.csv", index=False)
    print(f"\nПОПАДАНИЙ в тест: {len(m)}")
    print(m.to_string())


if __name__ == "__main__":
    main()
