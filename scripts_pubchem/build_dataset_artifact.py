# -*- coding: utf-8 -*-
"""
Сборка финального артефакта-датасета в формате задания:
smiles, IC50, CC50, SI, doi, source.

Источники: data/chembl_vacv.csv (с IC50/CC50/DOI), data/pubchem_vacv.csv,
DOI для PubChem-строк — из data/dataset_with_doi.csv (результат enrich_doi_v2.py).

Дедупликация по каноническому SMILES (медиана). Выход: data/dataset_final.csv

Запуск: python scripts/build_dataset_artifact.py
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from features import canonicalize  # noqa: E402

DATA = ROOT / "data"


def main():
    ch = pd.read_csv(DATA / "chembl_vacv.csv")
    pc = pd.read_csv(DATA / "pubchem_vacv.csv")
    cols = ["smiles", "ic50_nM", "cc50_nM", "si", "doi", "source"]
    al = pd.concat([ch[cols], pc[cols]], ignore_index=True)

    # DOI из обогащённой версии (там PubChem-строки получили DOI через ChEMBL-маппинг)
    # ключ — КАНОНИЧЕСКИЙ smiles, т.к. в dataset_with_doi.csv smiles уже канонизирован
    enr = pd.read_csv(DATA / "dataset_with_doi.csv")
    doi_map = {canonicalize(s): d for s, d in zip(enr["smiles"], enr["doi"])}

    al["canon"] = al["smiles"].map(canonicalize)
    al = al.dropna(subset=["canon"])
    al = al[al["si"] > 0]
    al["doi"] = [d if isinstance(d, str) and d.strip() else doi_map.get(c, "")
                 for d, c in zip(al["doi"], al["canon"])]

    out = (al.groupby("canon")
             .agg(ic50_nM=("ic50_nM", "median"),
                  cc50_nM=("cc50_nM", "median"),
                  si=("si", "median"),
                  doi=("doi", lambda s: next((x for x in s if isinstance(x, str) and x.strip()), "")),
                  source=("source", lambda s: ";".join(sorted(set(s)))),
                  n_measurements=("si", "size"))
             .reset_index()
             .rename(columns={"canon": "smiles"})
             .sort_values("si", ascending=False))
    out = out[["smiles", "ic50_nM", "cc50_nM", "si", "doi", "source", "n_measurements"]]
    out.to_csv(DATA / "dataset_final.csv", index=False)

    have_doi = (out["doi"].astype(str).str.strip() != "").sum()
    have_pair = out["ic50_nM"].notna().sum()
    print(f"Итоговый датасет: {len(out)} уникальных молекул -> data/dataset_final.csv")
    print(f"  с DOI: {have_doi} ({100*have_doi/len(out):.0f}%)")
    print(f"  с измеренными IC50+CC50: {have_pair}")
    print(f"  только с готовым SI (прямые измерения селективности): {len(out)-have_pair}")


if __name__ == "__main__":
    main()
