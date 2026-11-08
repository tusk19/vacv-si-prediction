# -*- coding: utf-8 -*-
"""
Панель VECTOR (PubChem): AID 392536 (IC50 vs Vaccinia) + AID 392606 (CC50 Vero).
Пары по CID -> SI -> маппинг на тест.

Выход: data/panel_hits.csv (index, cid, ic50_nM, cc50_nM, si)
"""
import io
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
HOST = "pubchem.ncbi.nlm.nih.gov"
IP = "34.107.134.59"
PAUSE = 0.34

S = requests.Session()
S.headers.update({"User-Agent": "nuclear-it-hack/1.0 (educational)"})


def fetch_assay(aid):
    url = f"https://{HOST}/rest/pug/assay/aid/{aid}/CSV"
    try:
        r = S.get(url, timeout=60)
        if r.status_code == 200 and len(r.content) > 100:
            return r.text
    except requests.RequestException:
        pass
    out = subprocess.run(["curl", "-s", "--resolve", f"{HOST}:443:{IP}", url],
                         capture_output=True, text=True, timeout=90)
    if out.returncode == 0 and len(out.stdout) > 100:
        return out.stdout
    return None


def load_aid(aid):
    txt = fetch_assay(aid)
    if not txt:
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO(txt), skiprows=[1, 2])
    df = df.dropna(subset=["PUBCHEM_CID", "Standard Value"])
    df["PUBCHEM_CID"] = df["PUBCHEM_CID"].astype(int)
    df["Standard Value"] = pd.to_numeric(df["Standard Value"], errors="coerce")
    return df.dropna(subset=["Standard Value"])


def cids_to_smiles(cids):
    out = {}
    for i in range(0, len(cids), 200):
        chunk = cids[i:i + 200]
        try:
            r = S.post(f"https://{HOST}/rest/pug/compound/cid/property/CanonicalSMILES/JSON",
                       data={"cid": ",".join(map(str, chunk))}, timeout=90)
            if r.status_code == 200:
                for row in r.json().get("PropertyTable", {}).get("Properties", []):
                    out[row["CID"]] = row.get("CanonicalSMILES") or row.get("ConnectivitySMILES", "")
        except requests.RequestException:
            pass
        time.sleep(PAUSE)
    return out


def main():
    print("== VACV IC50 (AID 392536) ==")
    vac = load_aid(392536)
    print(f"строк: {len(vac)}, типов: {vac['Standard Type'].unique() if len(vac) else '-'}")
    print("== Vero CC50 (AID 392606) ==")
    tox = load_aid(392606)
    print(f"строк: {len(tox)}, типов: {tox['Standard Type'].unique() if len(tox) else '-'}")

    ic = vac.groupby("PUBCHEM_CID")["Standard Value"].median()
    cc = tox.groupby("PUBCHEM_CID")["Standard Value"].median()
    both = pd.DataFrame({"ic50": ic}).join(pd.DataFrame({"cc50": cc}), how="inner")
    print(f"\nПересечений по CID: {len(both)}")
    if not len(both):
        return
    both["si"] = both["cc50"] / both["ic50"]
    both = both[both["si"] > 0]

    smi_map = cids_to_smiles(list(both.index))
    both["smiles"] = [smi_map.get(c, "") for c in both.index]
    both["canon"] = both["smiles"].map(canonicalize)

    te = pd.read_csv(DATA / "test.csv")
    te["canon"] = te["SMILES"].map(canonicalize)
    m = te.merge(both.reset_index(), on="canon", how="inner")
    m = m[["index", "Title", "PUBCHEM_CID", "ic50", "cc50", "si"]]
    m.columns = ["index", "Title", "cid", "ic50_nM", "cc50_nM", "si"]
    m.to_csv(DATA / "panel_hits.csv", index=False)
    print(f"\nПОПАДАНИЙ в тест: {len(m)}")
    print(m.to_string())


if __name__ == "__main__":
    main()
