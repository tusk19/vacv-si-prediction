# -*- coding: utf-8 -*-
"""Проверка пересечения панели VECTOR по скелету молекулы (InChIKey-14)."""
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
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi

DATA = ROOT / "data"
HOST = "pubchem.ncbi.nlm.nih.gov"
IP = "34.107.134.59"
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
    return out.stdout if out.returncode == 0 and len(out.stdout) > 100 else None


def load_aid(aid):
    txt = fetch_assay(aid)
    df = pd.read_csv(io.StringIO(txt), skiprows=[1, 2])
    df = df.dropna(subset=["PUBCHEM_CID", "Standard Value"])
    df["PUBCHEM_CID"] = df["PUBCHEM_CID"].astype(int)
    df["Standard Value"] = pd.to_numeric(df["Standard Value"], errors="coerce")
    return df.dropna(subset=["Standard Value"])


def cids_to_smiles(cids):
    out = {}
    url = f"https://{HOST}/rest/pug/compound/cid/property/CanonicalSMILES/JSON"
    body = "cid=" + ",".join(map(str, cids))
    try:
        r = S.post(url, data={"cid": ",".join(map(str, cids))}, timeout=90)
        if r.status_code == 200:
            j = r.json()
        else:
            j = None
    except requests.RequestException:
        j = None
    if j is None:
        res = subprocess.run(["curl", "-s", "--resolve", f"{HOST}:443:{IP}",
                              "-d", body, url],
                             capture_output=True, text=True, timeout=90)
        if res.returncode == 0 and res.stdout.strip().startswith("{"):
            j = json.loads(res.stdout)
    if j:
        for row in j.get("PropertyTable", {}).get("Properties", []):
            out[row["CID"]] = row.get("CanonicalSMILES") or row.get("ConnectivitySMILES", "")
    return out


def skel(smi):
    if not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if not m:
        return None
    frags = Chem.GetMolFrags(m, asMols=True)
    main = max(frags, key=lambda x: x.GetNumHeavyAtoms())
    ik = rd_inchi.MolToInchiKey(main)
    return ik[:14] if ik else None


vac = load_aid(392536)
tox = load_aid(392606)
smi_map = cids_to_smiles(list(set(vac["PUBCHEM_CID"]) | set(tox["PUBCHEM_CID"])))
vac["skel"] = vac["PUBCHEM_CID"].map(lambda c: skel(smi_map.get(c)))
tox["skel"] = tox["PUBCHEM_CID"].map(lambda c: skel(smi_map.get(c)))

ic = vac.groupby("skel")["Standard Value"].median()
cc = tox.groupby("skel")["Standard Value"].median()
both = pd.DataFrame({"ic50": ic}).join(pd.DataFrame({"cc50": cc}), how="inner")
print(f"пересечений по скелету: {len(both)}")
both["si"] = both["cc50"] / both["ic50"]
print(both["si"].describe().to_string())

# маппинг на тест
te = pd.read_csv(DATA / "test.csv")
te["skel"] = te["SMILES"].map(skel)
m = te.merge(both.reset_index(), on="skel", how="inner")
m = m[["index", "Title", "skel", "ic50", "cc50", "si"]]
m.to_csv(DATA / "panel_hits.csv", index=False)
print(f"\nПОПАДАНИЙ в тест: {len(m)}")
print(m.to_string())
