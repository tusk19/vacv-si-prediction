# -*- coding: utf-8 -*-
"""
Поиск покрытия по НАЗВАНИЯМ молекул (колонка Title в test.csv).

Логика: структурный матчинг мог промахнуться из-за другой соли/таутомера.
Ищем название в PubChem -> CID -> канонический SMILES -> проверяем,
есть ли такой SMILES в наших картах SI (paired_docs, combined, pubchem_map).
Заодно сохраняем CID для дальнейшего копания.

Вход:  data/uncovered_test.csv
Выход: data/name_search_results.csv (index, Title, cid, found_si, src)

Запуск: python scripts/search_by_name.py  (фоновая джоба, ~10-20 мин)
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from features import canonicalize  # noqa: E402

DATA = ROOT / "data"
PUG_HOST = "pubchem.ncbi.nlm.nih.gov"
PUG_IP = "34.107.134.59"

S = requests.Session()
S.headers.update({"User-Agent": "nuclear-it-hack/1.0 (educational)"})


def pug(url):
    try:
        r = S.get(url, timeout=30)
        if r.status_code == 200:
            return r.json()
    except requests.RequestException:
        pass
    try:
        out = subprocess.run(
            ["curl", "-s", "--resolve", f"{PUG_HOST}:443:{PUG_IP}", url],
            capture_output=True, text=True, timeout=45)
        if out.returncode == 0 and out.stdout.strip().startswith("{"):
            return json.loads(out.stdout)
    except Exception:
        pass
    return None


def name_to_smiles(name):
    q = requests.utils.quote(str(name).strip())
    d = pug(f"https://{PUG_HOST}/rest/pug/compound/name/{q}/property/CanonicalSMILES/JSON")
    if not d:
        return None, None
    try:
        p = d["PropertyTable"]["Properties"][0]
        return p.get("CID"), (p.get("CanonicalSMILES") or p.get("ConnectivitySMILES"))
    except (KeyError, IndexError):
        return None, None


def main():
    unc = pd.read_csv(DATA / "uncovered_test.csv")

    # карты SI из coverage_map-источников
    paired = pd.read_csv(DATA / "chembl_paired_docs.csv")
    paired = paired[paired["cc50_censored"] != True]
    paired["canon"] = paired["smiles"].map(canonicalize)
    m_paired = paired.groupby("canon")["si"].median().to_dict()
    comb = pd.read_csv(DATA / "combined_vacv.csv")
    comb["canon"] = comb["smiles"].map(canonicalize)
    m_comb = comb.groupby("canon")["si"].median().to_dict()
    pcm = pd.read_csv(DATA / "pubchem_si_map.csv")
    m_pc = pcm.groupby("canon")["si_pc"].median().to_dict()

    results = []
    out_path = DATA / "name_search_results.csv"
    for n, (_, row) in enumerate(unc.iterrows()):
        title = str(row["Title"]).strip()
        cid, smi = name_to_smiles(title)
        found_si, src = None, ""
        canon = canonicalize(smi) if smi else None
        if canon:
            if canon in m_paired:
                found_si, src = m_paired[canon], "paired_docs"
            elif canon in m_comb:
                found_si, src = m_comb[canon], "combined"
            elif canon in m_pc:
                found_si, src = m_pc[canon], "pubchem_map"
        results.append({"index": row["index"], "Title": title, "cid": cid,
                        "canon_via_name": canon, "found_si": found_si, "src": src})
        if n % 25 == 0:
            print(f"  {n}/{len(unc)} | найдено подстановок: "
                  f"{sum(1 for r in results if r['found_si'] is not None)}", flush=True)
            pd.DataFrame(results).to_csv(out_path, index=False)  # инкрементально
        time.sleep(0.35)

    df = pd.DataFrame(results)
    df.to_csv(out_path, index=False)
    hits = df[df["found_si"].notna()]
    print(f"\nГотово: {len(df)} названий проверено, CID найден для {df['cid'].notna().sum()},"
          f" подстановок SI: {len(hits)} -> {out_path}")


if __name__ == "__main__":
    main()
