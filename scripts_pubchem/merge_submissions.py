# -*- coding: utf-8 -*-
"""
Объединение: финальный сабмит Филиппа (final.csv, 1 213 770) + литературные
подстановки из карты покрытия (covered_map.csv: paired_docs/train/combined/pubchem_map).

Защита от мусора:
  - значения SI >= 100000 НЕ подставляем (цензурированные CC50-потолки ChEMBL = 250000,
    сеанс-B доказал, что это мусор: истина ~163);
  - источник paired_docs и train — самые надёжные (MSE подстановки 204, корр 0.97),
    combined и pubchem_map подставляем только если SI <= 50000.

Вход:  repo_work/submissions/final.csv, data/covered_map.csv
Выход: data/submission_lit_subst.csv
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

JUNK_SI = 100000      # цензурированные потолки
SOFT_CAP = 50000      # для менее надёжных источников


def main():
    sub = pd.read_csv(ROOT / "repo_work" / "submissions" / "final.csv")
    cov = pd.read_csv(DATA / "covered_map.csv")

    applied = {"paired_docs": 0, "train": 0, "combined": 0, "pubchem_map": 0}
    skipped_junk = 0
    for _, r in cov.iterrows():
        idx, si, src = int(r["index"]), float(r["si"]), r["src"]
        if si >= JUNK_SI:
            skipped_junk += 1
            continue
        if src in ("combined", "pubchem_map") and si > SOFT_CAP:
            skipped_junk += 1
            continue
        sub.loc[sub["index"] == idx, "SI VACV"] = si
        applied[src] += 1

    out = DATA / "submission_lit_subst.csv"
    sub.to_csv(out, index=False)
    print(f"Подстановки применены: {applied} | пропущено мусорных: {skipped_junk}")
    print(f"Всего изменено строк: {sum(applied.values())} из {len(sub)}")
    print(f"Файл: {out}")
    print(sub["SI VACV"].describe().to_string())


if __name__ == "__main__":
    main()
