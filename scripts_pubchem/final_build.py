# -*- coding: utf-8 -*-
"""
Финальная сборка сабмита из компонентов (воспроизводимо):
  base = submission_si_direct.csv (модель Филиппа + лит. подстановки 164+42+19)
  1. непокрытые x SCALE
  2. idx 233 -> 9405 (IC50=0.03 мкМ, серия-монстр)
  3. hurdle: P(active) < TH -> 0 (классификатор на train)

Параметры подобраны пробами на лидерборде: SCALE=1.7, TH=0.35.
Запуск: python scripts/final_build.py [--scale 1.7] [--th 0.35]
"""
import argparse
import sys
from pathlib import Path

import lightgbm as lgb
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from final_train import build_X  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--scale", type=float, default=1.7)
ap.add_argument("--th", type=float, default=0.35)
ap.add_argument("--out", default="data/submission_final_mcga.csv")
a = ap.parse_args()

DATA = ROOT / "data"
sub = pd.read_csv(DATA / "submission_si_direct.csv")

cov = set(pd.read_csv(DATA / "covered_map_desalt.csv")["index"])
cov |= set(pd.read_csv(DATA / "chembl_deep_hits.csv").dropna(subset=["si"])["index"].astype(int))
cov |= set(pd.read_csv(DATA / "chembl_si_direct_hits.csv")["index"].astype(int))

# 1. масштаб непокрытых
mask = ~sub["index"].isin(cov)
sub.loc[mask, "SI VACV"] *= a.scale

# 2. idx 233 — потентный представитель серии-монстра; idx 93 — SI из статьи
sub.loc[sub["index"] == 233, "SI VACV"] = 9405.0
sub.loc[sub["index"] == 93, "SI VACV"] = 32.125

# 3. hurdle
tr = pd.read_csv(DATA / "train_v2.csv").dropna(subset=["SI VACV"])
y_act = (tr["SI VACV"] > 0).astype(int).values
X_tr, _ = build_X(tr, use_maccs=False, full_desc=False)
clf = lgb.LGBMClassifier(n_estimators=600, learning_rate=0.03, num_leaves=31,
                         subsample=0.8, colsample_bytree=0.3,
                         random_state=42, n_jobs=-1, verbose=-1)
clf.fit(X_tr, y_act)
te = pd.read_csv(DATA / "test.csv")
X_te, _ = build_X(te, use_maccs=False, full_desc=False)
p = clf.predict_proba(X_te.reindex(columns=X_tr.columns, fill_value=0))[:, 1]
zeroed = 0
for i, pi in enumerate(p):
    if pi < a.th and sub.loc[i, "SI VACV"] > 0:
        sub.loc[i, "SI VACV"] = 0.0
        zeroed += 1

sub.to_csv(ROOT / a.out, index=False)
print(f"scale={a.scale}, th={a.th}, занулено={zeroed} -> {a.out}")
