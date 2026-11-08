# -*- coding: utf-8 -*-
"""
Проверка: есть ли в тесте молекулы, похожие (Tanimoto >= 0.7) на train-монстров
(SI > 2000). Для них, возможно, правильнее подставить SI монстра, а не потолок 5000.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

tr = pd.read_csv(DATA / "train_v2.csv").dropna(subset=["SI VACV"])
te = pd.read_csv(DATA / "test.csv")
sub = pd.read_csv(DATA / "submission_si_direct.csv")

monsters = tr[tr["SI VACV"] > 2000]
print(f"Монстров в train: {len(monsters)}")
print(monsters[["Title", "SI VACV"]].to_string())


def fp(s):
    m = Chem.MolFromSmiles(s)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None


mfps = [(fp(s), si, t) for s, si, t in zip(monsters["SMILES"], monsters["SI VACV"], monsters["Title"])]
mfps = [x for x in mfps if x[0] is not None]

sub_map = dict(zip(sub["index"], sub["SI VACV"]))
print("\nПохожие на монстров в тесте (Tanimoto >= 0.5):")
for _, r in te.iterrows():
    f = fp(r["SMILES"])
    if f is None:
        continue
    for mfp, si, title in mfps:
        sim = DataStructs.TanimotoSimilarity(f, mfp)
        if sim >= 0.5:
            idx = int(r["index"])
            print(f"  test idx {idx} ({r['Title']}): sim={sim:.2f} к train '{title}' "
                  f"(SI={si}) | наш прогноз: {sub_map.get(idx, '?'):.0f}")
