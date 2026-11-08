# -*- coding: utf-8 -*-
"""Валидация метода chembl_deep на train: прогоняем молекулы train через ту же
процедуру (similarity 100 -> VACV IC50 min + same-doc CC50 first) и сравниваем
подставленный SI с настоящим SI VACV.

Выход: MSE/MAE подстановки и по строкам — только для пересечений.
"""
import sys
import time
from pathlib import Path

import numpy as np
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


def deep_si(smiles):
    d = get(f"{BASE}/similarity/{requests.utils.quote(smiles, safe='')}/100.json")
    time.sleep(PAUSE)
    mols = (d or {}).get("molecules", [])
    if not mols:
        return None
    mid = mols[0]["molecule_chembl_id"]
    acts = get(f"{BASE}/activity.json",
               {"molecule_chembl_id": mid, "target_chembl_id": VACV, "limit": 100})
    time.sleep(PAUSE)
    ics = [a for a in (acts or {}).get("activities", [])
           if a.get("standard_type") in ("IC50", "EC50")
           and a.get("standard_value") and a.get("standard_units") == "nM"]
    if not ics:
        return None
    a0 = min(ics, key=lambda a: float(a["standard_value"]))
    ic_uM = float(a0["standard_value"]) / 1000.0
    doc = a0.get("document_chembl_id")
    if not doc:
        return None
    d2 = get(f"{BASE}/activity.json",
             {"document_chembl_id": doc, "molecule_chembl_id": mid,
              "standard_type": "CC50", "limit": 100})
    time.sleep(PAUSE)
    ccs = [a for a in (d2 or {}).get("activities", [])
           if a.get("standard_value") and a.get("standard_units") == "nM"]
    if not ccs or ic_uM <= 0:
        return None
    return float(ccs[0]["standard_value"]) / 1000.0 / ic_uM


def main():
    tr = pd.read_csv(DATA / "train_v2.csv").dropna(subset=["SI VACV"])
    tr = tr[tr["SI VACV"] > 0].head(60)  # достаточно для оценки
    rows = []
    for n, (_, r) in enumerate(tr.iterrows()):
        si = deep_si(r["SMILES"])
        if si:
            rows.append({"true": r["SI VACV"], "subst": si})
        if (n + 1) % 10 == 0:
            print(f"  {n+1}/{len(tr)} | покрыто {len(rows)}", flush=True)
    df = pd.DataFrame(rows)
    if not len(df):
        print("пересечений нет")
        return
    mse = float(np.mean((df["true"] - df["subst"]) ** 2))
    mae = float(np.mean(abs(df["true"] - df["subst"])))
    ratio = (df["subst"] / df["true"]).replace([np.inf], np.nan).dropna()
    print(f"\nПересечений: {len(df)}")
    print(f"MSE подстановки: {mse:.1f} | MAE: {mae:.2f}")
    print(f"медиана отношения subst/true: {ratio.median():.3f}")
    print(f"доля с ошибкой <2x: {((ratio > 0.5) & (ratio < 2)).mean():.0%}")


if __name__ == "__main__":
    main()
