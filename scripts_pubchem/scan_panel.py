# -*- coding: utf-8 -*-
"""Скан диапазона PubChem AID вокруг 3925xx: ищем тесты vaccinia/orthopox и CC50/cytotox."""
import json
import subprocess
import time

HOST = "pubchem.ncbi.nlm.nih.gov"
IP = "34.107.134.59"


def desc(aid):
    url = f"https://{HOST}/rest/pug/assay/aid/{aid}/description/JSON"
    try:
        out = subprocess.run(["curl", "-s", "--resolve", f"{HOST}:443:{IP}", url],
                             capture_output=True, text=True, timeout=45)
        if out.returncode == 0 and out.stdout.strip().startswith("{"):
            d = json.loads(out.stdout)
            a = d["PC_AssayContainer"][0]["assay"]["descr"]
            return a.get("name", "")
    except Exception:
        pass
    return None


found = []
for aid in range(392495, 392580):
    name = desc(aid)
    if name:
        low = name.lower()
        if any(k in low for k in ("vaccinia", "orthopox", "cytotox", "cc50", "toxicity", "variola", "pox")):
            print(aid, "|", name[:120])
            found.append(aid)
    time.sleep(0.3)
print("\nнайдено:", found)
