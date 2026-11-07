# -*- coding: utf-8 -*-
"""
Скачивание данных об активности молекул против вируса оспы (VACV) из ChEMBL.

Что делает:
1. Тянет ВСЕ записи об активности против target CHEMBL613493 (Vaccinia virus).
2. Для каждой молекулы ищет CC50 (цитотоксичность) по всему ChEMBL.
3. Считает SI = CC50 / IC50.
4. Сохраняет CSV: smiles, ic50, cc50, si, doi, source.

Запуск:  python scripts/download_chembl.py
Выход:   data/chembl_vacv.csv  +  data/chembl_vacv_raw.csv (все сырые записи)
"""
import json
import time
import sys
from pathlib import Path
from urllib.parse import urlencode

import requests
import pandas as pd

BASE = "https://www.ebi.ac.uk/chembl/api/data"
VACV_TARGET = "CHEMBL613493"  # Vaccinia virus (organism)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "nuclear-it-hack/1.0 (educational)"})


def fetch_all(url, params=None, pause=0.4, label=""):
    """Выкачивает все страницы пагинированного ответа ChEMBL API."""
    out = []
    while url:
        for attempt in range(4):
            try:
                r = SESSION.get(url, params=params, timeout=60)
                r.raise_for_status()
                break
            except requests.RequestException as e:
                print(f"  retry {attempt+1}: {e}", flush=True)
                time.sleep(2 * (attempt + 1))
        else:
            raise RuntimeError(f"не удалось скачать {url}")
        d = r.json()
        key = [k for k in d if k != "page_meta"][0]
        out.extend(d[key])
        nxt = d["page_meta"].get("next")
        if nxt:
            # next приходит как путь вида /chembl/api/data/... — доклеиваем только хост
            url = nxt if nxt.startswith("http") else ("https://www.ebi.ac.uk" + nxt)
        else:
            url = None
        params = None  # next-ссылка уже содержит параметры
        total = d["page_meta"].get("total_count", "?")
        print(f"  {label}: скачано {len(out)}/{total}", flush=True)
        time.sleep(pause)
    return out


def main():
    # ---------- 1. Активность против VACV ----------
    print("== Шаг 1. Записи активности против Vaccinia virus ==", flush=True)
    acts = fetch_all(
        f"{BASE}/activity.json",
        params={"target_chembl_id": VACV_TARGET, "limit": 1000},
        label="VACV",
    )
    print(f"Всего записей: {len(acts)}")

    rows = []
    for a in acts:
        rows.append({
            "molecule_chembl_id": a.get("molecule_chembl_id"),
            "smiles": (a.get("canonical_smiles") or "").strip(),
            "standard_type": a.get("standard_type"),
            "standard_value": a.get("standard_value"),
            "standard_units": a.get("standard_units"),
            "standard_relation": a.get("standard_relation"),
            "assay_type": a.get("assay_type"),
            "assay_description": (a.get("assay_description") or "")[:200],
            "document_chembl_id": a.get("document_chembl_id"),
        })
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "chembl_vacv_raw.csv", index=False)
    print("\nРаспределение типов измерений против VACV:")
    print(df["standard_type"].value_counts().head(15).to_string())

    # ---------- 2. DOI документов ----------
    print("\n== Шаг 2. DOI статей ==", flush=True)
    doc_ids = sorted({d for d in df["document_chembl_id"] if d})
    doi_map = {}
    # ChEMBL API позволяет фильтровать списком id через __in (пачками)
    for i in range(0, len(doc_ids), 200):
        chunk = doc_ids[i:i + 200]
        docs = fetch_all(
            f"{BASE}/document.json",
            params={"document_chembl_id__in": ",".join(chunk), "limit": 1000},
            label="docs",
        )
        for d in docs:
            doi_map[d["document_chembl_id"]] = d.get("doi") or ""
    print(f"DOI найдено для {sum(1 for v in doi_map.values() if v)} из {len(doc_ids)} документов")

    # ---------- 3. IC50 против вируса ----------
    ic50 = df[(df["standard_type"].isin(["IC50", "EC50"])) &
              (df["standard_units"] == "nM") &
              (df["standard_relation"].isin(["=", "<", "<="])) &
              (df["smiles"] != "")].copy()
    ic50["standard_value"] = pd.to_numeric(ic50["standard_value"], errors="coerce")
    ic50 = ic50.dropna(subset=["standard_value"])
    # усредняем дубликаты по молекуле (медиана устойчивее к выбросам)
    ic50_agg = (ic50.groupby("smiles")
                    .agg(ic50_nM=("standard_value", "median"),
                         n_ic50=("standard_value", "size"),
                         doc=("document_chembl_id", "first"))
                    .reset_index())
    print(f"\nМолекул с IC50/EC50 против VACV: {len(ic50_agg)}")

    # ---------- 4. CC50 (цитотоксичность) для тех же молекул ----------
    print("\n== Шаг 3. Поиск CC50 для этих молекул по всему ChEMBL ==", flush=True)
    mol_ids = sorted(set(df[df["smiles"].isin(ic50_agg["smiles"])]["molecule_chembl_id"]))
    cc50_rows = []
    for i in range(0, len(mol_ids), 50):
        chunk = mol_ids[i:i + 50]
        res = fetch_all(
            f"{BASE}/activity.json",
            params={"molecule_chembl_id__in": ",".join(chunk),
                    "standard_type": "CC50", "limit": 1000},
            label=f"CC50 {i}-{i+len(chunk)}",
        )
        for a in res:
            if (a.get("standard_units") == "nM" and a.get("canonical_smiles")
                    and a.get("standard_value") is not None
                    and a.get("standard_relation") in ("=", ">", ">=")):
                cc50_rows.append({
                    "smiles": a["canonical_smiles"],
                    "cc50_value": float(a["standard_value"]),
                })
    cc50_df = pd.DataFrame(cc50_rows)
    if len(cc50_df):
        cc50_agg = (cc50_df.groupby("smiles")
                           .agg(cc50_nM=("cc50_value", "median"),
                                n_cc50=("cc50_value", "size"))
                           .reset_index())
    else:
        cc50_agg = pd.DataFrame(columns=["smiles", "cc50_nM", "n_cc50"])
    print(f"Молекул с найденным CC50: {len(cc50_agg)}")

    # ---------- 5. Склейка: SI = CC50 / IC50 ----------
    out = ic50_agg.merge(cc50_agg, on="smiles", how="inner")
    out["si"] = out["cc50_nM"] / out["ic50_nM"]
    out["doi"] = out["doc"].map(doi_map).fillna("")
    out["source"] = "ChEMBL:" + out["doc"].fillna("")
    out = out[["smiles", "ic50_nM", "cc50_nM", "si", "doi", "source", "n_ic50", "n_cc50"]]
    out = out.sort_values("si", ascending=False)
    out.to_csv(DATA / "chembl_vacv.csv", index=False)
    print(f"\n== Готово: {len(out)} молекул с парой IC50+CC50 -> data/chembl_vacv.csv ==")
    print(out["si"].describe().to_string())


if __name__ == "__main__":
    main()
