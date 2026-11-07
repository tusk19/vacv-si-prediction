# -*- coding: utf-8 -*-
"""
Обогащение датасета DOI (v2).

Путь: PubChem AID -> ChEMBL assay id (в описании теста лежит aid_source)
      -> ChEMBL document -> DOI.
Строки, напрямую пришедшие из ChEMBL, DOI уже имеют — не трогаем.

DNS pubchem периодически падает у провайдера: для PUG-запросов есть обход
через curl --resolve с IP из DoH.

Вход:  data/combined_vacv.csv
Выход: data/dataset_with_doi.csv  (smiles, ic50_nM, cc50_nM, si, doi, source)

Запуск: python scripts/enrich_doi_v2.py
"""
import json
import re
import subprocess
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PUG_IP = "34.107.134.59"  # pubchem.ncbi.nlm.nih.gov (IP из DoH cloudflare)
PUG_HOST = "pubchem.ncbi.nlm.nih.gov"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

S = requests.Session()
S.headers.update({"User-Agent": "nuclear-it-hack/1.0 (educational)"})


def pug_description(aid):
    """AID -> JSON описания теста. Сначала обычный запрос, при DNS-сбое — через IP."""
    url = f"https://{PUG_HOST}/rest/pug/assay/aid/{aid}/description/JSON"
    try:
        r = S.get(url, timeout=40)
        if r.status_code == 200:
            return r.json()
    except requests.RequestException:
        pass
    try:
        out = subprocess.run(
            ["curl", "-s", "--resolve", f"{PUG_HOST}:443:{PUG_IP}", url],
            capture_output=True, text=True, timeout=60)
        if out.returncode == 0 and out.stdout.strip().startswith("{"):
            return json.loads(out.stdout)
    except Exception:
        pass
    return None


def aid_to_chembl_assay(aid):
    d = pug_description(aid)
    if not d:
        return None
    try:
        descr = d["PC_AssayContainer"][0]["assay"]["descr"]
        src = descr.get("aid_source", {}).get("db", {})
        if src.get("name") == "ChEMBL":
            return src["source_id"]["str"]
    except (KeyError, IndexError, TypeError):
        pass
    return None


def chembl_assay_to_doi(assay_ids):
    """Пачка ChEMBL assay id -> {assay_id: doi}."""
    out = {}
    ids = sorted(set(assay_ids))
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        for attempt in range(4):
            try:
                r = S.get(f"{CHEMBL}/assay.json",
                          params={"assay_chembl_id__in": ",".join(chunk), "limit": 1000},
                          timeout=60)
                r.raise_for_status()
                break
            except requests.RequestException:
                time.sleep(2 * (attempt + 1))
        else:
            continue
        doc_ids = {}
        for a in r.json().get("assays", []):
            if a.get("document_chembl_id"):
                doc_ids[a["assay_chembl_id"]] = a["document_chembl_id"]
        docs = sorted(set(doc_ids.values()))
        if docs:
            for attempt in range(4):
                try:
                    r2 = S.get(f"{CHEMBL}/document.json",
                               params={"document_chembl_id__in": ",".join(docs), "limit": 1000},
                               timeout=60)
                    r2.raise_for_status()
                    break
                except requests.RequestException:
                    time.sleep(2 * (attempt + 1))
            else:
                continue
            d2doi = {d["document_chembl_id"]: (d.get("doi") or "")
                     for d in r2.json().get("documents", [])}
            for a_id, did in doc_ids.items():
                out[a_id] = d2doi.get(did, "")
        print(f"  ChEMBL: {min(i + 50, len(ids))}/{len(ids)}", flush=True)
        time.sleep(0.4)
    return out


def main():
    df = pd.read_csv(DATA / "combined_vacv.csv")
    print(f"Строк: {len(df)}")

    def extract_aid(src):
        m = re.search(r"AID(\d+)", str(src))
        return m.group(1) if m else None

    df["aid"] = df["source"].map(extract_aid)
    aids = sorted({a for a in df["aid"] if isinstance(a, str) and a})
    print(f"Уникальных PubChem AID: {len(aids)}")

    print("== AID -> ChEMBL assay ==", flush=True)
    a2c = {}
    for n, aid in enumerate(aids):
        c = aid_to_chembl_assay(aid)
        if c:
            a2c[aid] = c
        if n % 10 == 0:
            print(f"  {n}/{len(aids)}", flush=True)
        time.sleep(0.3)
    print(f"Маппингов найдено: {len(a2c)}")

    print("== ChEMBL assay -> DOI ==", flush=True)
    c2doi = chembl_assay_to_doi(list(a2c.values()))

    def fill_doi(row):
        if isinstance(row["doi"], str) and row["doi"].strip():
            return row["doi"]
        if row["aid"] and row["aid"] in a2c:
            return c2doi.get(a2c[row["aid"]], "")
        return ""

    df["doi"] = df.apply(fill_doi, axis=1)
    out = df.drop(columns=["aid"])
    out.to_csv(DATA / "dataset_with_doi.csv", index=False)
    have = (out["doi"].astype(str).str.strip() != "").sum()
    print(f"\nГотово: {len(out)} строк -> data/dataset_with_doi.csv")
    print(f"Строк с DOI: {have} ({100 * have / len(out):.0f}%)")


if __name__ == "__main__":
    main()
