# -*- coding: utf-8 -*-
"""Добавляем глубокие покрытия (ChEMBL deep + PubChem deep) к лучшему сабмиту.

Вход:  data/submission_lit_subst.csv (текущий лучший, 1 182 246)
       data/chembl_deep_hits.csv, data/deep_coverage_hits.csv (если есть)
Выход: data/submission_deep_subst.csv
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
JUNK_SI = 100000


def main():
    sub = pd.read_csv(DATA / "submission_lit_subst.csv")
    applied = 0

    for f, src in [("chembl_deep_hits.csv", "chembl_deep"),
                   ("deep_coverage_hits.csv", "pubchem_deep")]:
        p = DATA / f
        if not p.exists():
            print(f"{f} — нет файла, пропуск")
            continue
        try:
            hits = pd.read_csv(p).dropna(subset=["si"])
        except pd.errors.EmptyDataError:
            print(f"{f} — пусто, пропуск")
            continue
        hits = hits[(hits["si"] > 0) & (hits["si"] < JUNK_SI)]
        # медиана по молекуле на случай нескольких находок
        hits = hits.groupby("index")["si"].median()
        for idx, si in hits.items():
            sub.loc[sub["index"] == int(idx), "SI VACV"] = si
            applied += 1
        print(f"{src}: подставлено {len(hits)}")

    out = DATA / "submission_deep_subst.csv"
    sub.to_csv(out, index=False)
    print(f"\nВсего новых подстановок: {applied}")
    print(f"Файл: {out}")
    print(sub["SI VACV"].describe().to_string())


if __name__ == "__main__":
    main()
