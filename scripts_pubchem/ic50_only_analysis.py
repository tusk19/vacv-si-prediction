# -*- coding: utf-8 -*-
"""
Анализ ic50_only молекул: тянем их IC50 против VACV и ищем CC50 в PubChem.
Потенциальные монстры = очень малый IC50 (нМ) без известного CC50.

Выход: data/ic50_only_enriched.csv (index, chembl_id, ic50_nM, cc50_nM_pubchem, si)
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
BASE = "https://www.ebi.ac.uk/chembl/api/data"
VACV = "CHEMBL613493"
PUG_HOST = "pubchem.ncbi.nlm.nih.gov"
PUG_IP = "34.107.134.59"
PAUSE = 0.34

S = requests.Session()
S.headers.update({"User-Agent": "nuclear-it-hack/1.0 (educational)"})


def get(url, params=None):
    for attempt in range(4):
        try:
            r = S.get(url, params=params, timeout=60)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
        time.sleep(1.5 * (attempt + 1))
    return None


def pug(url):
    try:
        r = S.get(url, timeout=40)
        if r.status_code == 200:
            return r.json()
    except requests.RequestException:
        pass
    try:
        out = subprocess.run(["curl", "-s", "--resolve", f"{PUG_HOST}:443:{PUG_IP}", url],
                             capture_output=True, text=True, timeout=60)
        if out.returncode == 0 and out.stdout.strip().startswith("{"):
            return json.loads(out.stdout)
    except Exception:
        pass
    return None


def main():
    hits = pd.read_csv(DATA / "chembl_deep_hits.csv")
    only = hits[hits["kind"] == "ic50_only"].copy()
    print(f"ic50_only: {len(only)}")

    # 1. IC50 против VACV (минимальный) для каждой
    rows = []
    for _, row in only.iterrows():
        mid = row["chembl_id"]
        a = get(f"{BASE}/activity.json",
                {"molecule_chembl_id": mid, "target_chembl_id": VACV, "limit": 100})
        time.sleep(PAUSE)
        ics = [x for x in (a or {}).get("activities", [])
               if x.get("standard_type") in ("IC50", "EC50")
               and x.get("standard_value") and x.get("standard_units") == "nM"]
        ic50 = min((float(x["standard_value"]) for x in ics), default=None)
        rows.append({"index": row["index"], "chembl_id": mid, "ic50_nM": ic50})
    df = pd.DataFrame(rows).dropna(subset=["ic50_nM"])
    df = df.sort_values("ic50_nM")
    print(df.head(15).to_string())

    # 2. Для самых потентных (IC50 < 500 нМ) ищем CC50 в PubChem
    potent = df[df["ic50_nM"] < 500]
    print(f"\nПотентных (IC50<500 нМ): {len(potent)} — ищем CC50 в PubChem")
    out = []
    for _, row in potent.iterrows():
        mid = row["chembl_id"]
        # chembl -> smiles
        m = get(f"{BASE}/molecule/{mid}.json")
        time.sleep(PAUSE)
        smi = ((m or {}).get("molecule_structures") or {}).get("canonical_smiles")
        if not smi:
            continue
        q = requests.utils.quote(smi, safe="")
        d = pug(f"https://{PUG_HOST}/rest/pug/compound/smiles/{q}/cids/JSON")
        time.sleep(PAUSE)
        cids = (d or {}).get("IdentifierList", {}).get("CID", [])
        if not cids:
            continue
        cid = cids[0]
        # assays CID -> ищем CC50
        d2 = pug(f"https://{PUG_HOST}/rest/pug/compound/cid/{cid}/assaysummary/JSON")
        time.sleep(PAUSE)
        aids = []
        for r in (d2 or {}).get("Table", {}).get("Row", []):
            c = r.get("Cell", [])
            if len(c) >= 10 and any(k in str(c[9]).lower() for k in
                                    ("cytotox", "cc50", "vero", "cell viability", "toxicity")):
                aids.append(str(c[0]))
        cc50 = None
        for aid in aids[:4]:
            try:
                r = S.get(f"https://{PUG_HOST}/rest/pug/assay/aid/{aid}/CSV", timeout=60)
                if r.status_code != 200:
                    continue
                adf = pd.read_csv(io.StringIO(r.text), skiprows=[1, 2])
                adf = adf[adf.get("PUBCHEM_CID") == cid]
                st = adf.get("Standard Type", pd.Series(dtype=str)).astype(str).str.upper()
                vals = adf[st == "CC50"]["Standard Value"]
                vals = pd.to_numeric(vals, errors="coerce").dropna()
                if len(vals):
                    cc50 = float(vals.median())
                    break
            except Exception:
                continue
            finally:
                time.sleep(PAUSE)
        si = cc50 / row["ic50_nM"] if cc50 and row["ic50_nM"] else None
        out.append({"index": row["index"], "chembl_id": mid,
                    "ic50_nM": row["ic50_nM"], "cc50_nM_pubchem": cc50, "si": si})
        print(f"  idx {row['index']}: IC50={row['ic50_nM']:.1f} нМ, "
              f"CC50_pubchem={cc50}, SI={si}", flush=True)

    pd.DataFrame(out).to_csv(DATA / "ic50_only_enriched.csv", index=False)
    print(f"\nГотово -> data/ic50_only_enriched.csv")


if __name__ == "__main__":
    main()
