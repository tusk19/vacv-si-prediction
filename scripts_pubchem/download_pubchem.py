# -*- coding: utf-8 -*-
"""
Скачивание данных по Vaccinia virus (VACV) из PubChem BioAssay.

Стратегия:
1. Находим тесты (assays), которые НАПРЯМУЮ измеряют индекс селективности
   (запрос: "vaccinia virus" AND "selectivity index") — это готовые SI.
2. Отдельно тянем тесты с IC50/EC50 против вируса и CC50 на клетках
   (запрос: "vaccinia virus" AND (IC50 OR CC50)) и склеиваем по CID.
3. CID -> канонический SMILES через PUG-REST.
4. Выход: data/pubchem_vacv.csv  (smiles, ic50_nM, cc50_nM, si, doi, source)

Запуск:  python scripts/download_pubchem.py
"""
import io
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PAUSE = 0.35  # PubChem просит <=5 запросов/сек

S = requests.Session()
S.headers.update({"User-Agent": "nuclear-it-hack/1.0 (educational)"})


def get(url, **kw):
    for attempt in range(5):
        try:
            r = S.get(url, timeout=90, **kw)
            if r.status_code == 200:
                return r
            if r.status_code in (404, 400):
                return None
        except requests.RequestException:
            pass
        time.sleep(1.5 * (attempt + 1))
    return None


def esearch(term, retmax=1000):
    r = get(f"{EUTILS}/esearch.fcgi", params={
        "db": "pcassay", "term": term, "retmode": "json", "retmax": retmax})
    if not r:
        return []
    return r.json()["esearchresult"]["idlist"]


def fetch_assay(aid):
    """Скачивает таблицу результатов теста. Возвращает DataFrame или None."""
    r = get(f"{PUG}/assay/aid/{aid}/CSV")
    if r is None or len(r.content) < 100:
        return None
    try:
        df = pd.read_csv(io.StringIO(r.text), skiprows=[1, 2])  # строки 2-3 — описания колонок
    except Exception:
        return None
    if "PUBCHEM_CID" not in df.columns or "Standard Value" not in df.columns:
        return None
    df = df.dropna(subset=["PUBCHEM_CID", "Standard Value"])
    df["PUBCHEM_CID"] = df["PUBCHEM_CID"].astype(int)
    df["Standard Value"] = pd.to_numeric(df["Standard Value"], errors="coerce")
    df = df.dropna(subset=["Standard Value"])
    df["aid"] = aid
    return df


def cids_to_smiles(cids):
    """Пачками переводит PubChem CID -> канонический SMILES."""
    out = {}
    cids = sorted(set(cids))
    for i in range(0, len(cids), 200):
        chunk = cids[i:i + 200]
        r = S.post(f"{PUG}/compound/cid/property/CanonicalSMILES/JSON",
                   data={"cid": ",".join(map(str, chunk))}, timeout=90)
        if r is not None and r.status_code == 200:
            for row in r.json().get("PropertyTable", {}).get("Properties", []):
                out[row["CID"]] = row.get("CanonicalSMILES") or row.get("ConnectivitySMILES", "")
        print(f"  SMILES: {min(i + 200, len(cids))}/{len(cids)}", flush=True)
        time.sleep(PAUSE)
    return out


def main():
    # ---------- 1. Тесты с готовым SI ----------
    print("== Поиск тестов 'selectivity index' ==", flush=True)
    si_aids = esearch("vaccinia virus AND selectivity index")
    print(f"Найдено тестов: {len(si_aids)}")

    si_rows = []
    for n, aid in enumerate(si_aids):
        df = fetch_assay(aid)
        if df is not None and len(df):
            si_rows.append(df[["PUBCHEM_CID", "Standard Value", "aid"]])
        if n % 20 == 0:
            print(f"  SI-тесты: {n}/{len(si_aids)}", flush=True)
        time.sleep(PAUSE)

    si_df = pd.concat(si_rows, ignore_index=True) if si_rows else pd.DataFrame()
    print(f"Строк с готовым SI: {len(si_df)}")

    # ---------- 2. Тесты с IC50/CC50 (склеиваем по CID) ----------
    print("\n== Поиск тестов IC50/CC50 ==", flush=True)
    pot_aids = esearch("vaccinia virus AND (IC50 OR CC50)")
    pot_aids = [a for a in pot_aids if a not in si_aids]
    print(f"Найдено тестов: {len(pot_aids)}")

    pot_rows = []
    for n, aid in enumerate(pot_aids):
        df = fetch_assay(aid)
        if df is not None and len(df):
            st = df["Standard Type"].astype(str).str.upper().str.strip() if "Standard Type" in df else pd.Series("", index=df.index)
            for stype, frame in df.groupby(st):
                if stype in ("IC50", "EC50", "CC50"):
                    pot_rows.append(pd.DataFrame({
                        "cid": frame["PUBCHEM_CID"], "type": stype,
                        "value": frame["Standard Value"], "aid": aid}))
        if n % 20 == 0:
            print(f"  IC50/CC50-тесты: {n}/{len(pot_aids)}", flush=True)
        time.sleep(PAUSE)

    pot_df = pd.concat(pot_rows, ignore_index=True) if pot_rows else pd.DataFrame()
    print(f"Строк IC50/CC50: {len(pot_df)}")

    # ---------- 3. SMILES ----------
    print("\n== CID -> SMILES ==", flush=True)
    all_cids = set()
    if len(si_df):
        all_cids |= set(si_df["PUBCHEM_CID"])
    if len(pot_df):
        all_cids |= set(pot_df["cid"])
    smi_map = cids_to_smiles(list(all_cids))
    print(f"SMILES получено для {len(smi_map)} CID")

    # ---------- 4. Сборка результата ----------
    results = []

    # прямые SI
    if len(si_df):
        si_df["smiles"] = si_df["PUBCHEM_CID"].map(smi_map)
        g = (si_df.dropna(subset=["smiles"])
                  .groupby("smiles")
                  .agg(si=("Standard Value", "median"),
                       n=("Standard Value", "size"),
                       aid=("aid", "first"))
                  .reset_index())
        for _, r in g.iterrows():
            results.append({"smiles": r["smiles"], "ic50_nM": None, "cc50_nM": None,
                            "si": r["si"], "doi": "",
                            "source": f"PubChem:AID{r['aid']}"})

    # пары IC50+CC50 по CID
    if len(pot_df):
        pot_df["smiles"] = pot_df["cid"].map(smi_map)
        pv = (pot_df.dropna(subset=["smiles"])
                    .pivot_table(index="smiles", columns="type",
                                 values="value", aggfunc="median"))
        have_ic = [c for c in ("IC50", "EC50") if c in pv.columns]
        if have_ic and "CC50" in pv.columns:
            pv["ic50"] = pv[have_ic].median(axis=1)
            pair = pv.dropna(subset=["ic50", "CC50"]).reset_index()
            for _, r in pair.iterrows():
                if r["ic50"] > 0:
                    results.append({"smiles": r["smiles"], "ic50_nM": r["ic50"],
                                    "cc50_nM": r["CC50"], "si": r["CC50"] / r["ic50"],
                                    "doi": "", "source": "PubChem:paired"})

    out = pd.DataFrame(results)
    out = out[out["si"] > 0]
    out = out.sort_values("si", ascending=False)
    out.to_csv(DATA / "pubchem_vacv.csv", index=False)
    print(f"\n== Готово: {len(out)} молекул -> data/pubchem_vacv.csv ==")
    print(out["si"].describe().to_string())


if __name__ == "__main__":
    main()
