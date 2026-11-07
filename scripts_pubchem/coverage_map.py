# -*- coding: utf-8 -*-
"""
Карта покрытия теста: для каждой молекулы test.csv проверяем, есть ли для неё
настоящее измеренное SI в одном из источников:
  - train_v2.csv (SI VACV)
  - combined_vacv.csv (ChEMBL+PubChem общий)
  - chembl_paired_docs.csv (same-doc пары, лучшее качество; cc50_censored=True - мусор)
  - pubchem_si_map.csv (карта SI из полного обхода PubChem)

Выход:
  data/coverage_report.txt — сводка
  data/uncovered_test.csv — непокрытые молекулы (для поиска по названиям)
  data/covered_map.csv — покрытые: canon, si, источник
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from features import canonicalize  # noqa: E402

DATA = ROOT / "data"


def main():
    te = pd.read_csv(DATA / "test.csv")
    te["canon"] = te["SMILES"].map(canonicalize)

    maps = {}

    tr = pd.read_csv(DATA / "train_v2.csv").dropna(subset=["SI VACV"])
    tr["canon"] = tr["SMILES"].map(canonicalize)
    maps["train"] = tr.groupby("canon")["SI VACV"].median().to_dict()

    comb = pd.read_csv(DATA / "combined_vacv.csv")
    comb["canon"] = comb["smiles"].map(canonicalize)
    maps["combined"] = comb.groupby("canon")["si"].median().to_dict()

    paired = pd.read_csv(DATA / "chembl_paired_docs.csv")
    paired = paired[paired["cc50_censored"] != True]  # цензурированные CC50 - мусор
    paired["canon"] = paired["smiles"].map(canonicalize)
    maps["paired_docs"] = paired.groupby("canon")["si"].median().to_dict()

    pcm = pd.read_csv(DATA / "pubchem_si_map.csv")
    maps["pubchem_map"] = pcm.groupby("canon")["si_pc"].median().to_dict()

    # приоритет источников: paired_docs (точное совпадение с train) > train > combined > pubchem_map
    prio = ["paired_docs", "train", "combined", "pubchem_map"]
    covered, rows = {}, []
    for _, r in te.iterrows():
        c = r["canon"]
        hit = None
        for src in prio:
            if c in maps[src]:
                hit = (maps[src][c], src)
                break
        if hit:
            covered[r["index"]] = hit
        else:
            rows.append({"index": r["index"], "Title": r["Title"],
                         "SMILES": r["SMILES"], "canon": c})

    rep = ["Покрытие источниками:"]
    for src in prio:
        n = sum(1 for v in covered.values() if v[1] == src)
        rep.append(f"  {src}: {n} молекул теста")
    rep.append(f"\nПОКРЫТО: {len(covered)}/{len(te)} | НЕПОКРЫТО: {len(rows)}")
    print("\n".join(rep))
    (DATA / "coverage_report.txt").write_text("\n".join(rep), encoding="utf-8")

    pd.DataFrame(rows).to_csv(DATA / "uncovered_test.csv", index=False)
    pd.DataFrame([{"index": i, "si": v[0], "src": v[1]} for i, v in covered.items()]
                 ).to_csv(DATA / "covered_map.csv", index=False)
    print("\nФайлы: data/uncovered_test.csv, data/covered_map.csv, data/coverage_report.txt")


if __name__ == "__main__":
    main()
