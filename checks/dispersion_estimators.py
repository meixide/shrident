# Recompute pairwise vs DerSimonian-Laird dispersion for NSW and Card-Krueger
# from the full-run CSVs in data/. Run from the repo root:
#   python checks/dispersion_estimators.py
import numpy as np, pandas as pd

def dl(p, v):
    w = 1 / v; pb = (w * p).sum() / w.sum(); Q = ((p - pb) ** 2 / v).sum()
    C = w.sum() - (w ** 2).sum() / w.sum(); return max(0, (Q - (len(p) - 1)) / C)

def pairwise(p, v):
    i, j = np.triu_indices(len(p), 1); return max(0, np.mean(((p[i] - p[j]) ** 2 - (v[i] + v[j])) / 2))

def report(tag, p, v):
    out = []
    for f in (pairwise, dl):
        t = f(p, v); w = 1 / (v + t)
        out.append(f"{f.__name__}: tau={np.sqrt(t):.3f} psi_EB={(w * p).sum() / w.sum():.3f}")
    print(f"{tag:28s} J={len(p):2d} | " + " | ".join(out))

ck = pd.read_csv("data/ck94_functionals.csv")
for m in ["DIM", "OM", "IPW"]:
    s = ck[ck.method == m]; report(f"CK {m} only", s.psi.to_numpy(), s.v.to_numpy())
report("CK all", ck.psi.to_numpy(), ck.v.to_numpy())
f = pd.read_csv("data/lalonde_functionals_hapc.csv"); p, v = f.psi.to_numpy(), f.v.to_numpy()
adj = ~f.label.str.match(r"^DIM-(CPS|PSID)")
report("NSW trial only", p[:3], v[:3]); report("NSW covariate-adjusted", p[adj], v[adj]); report("NSW all", p, v)
# Working-independence SE reported in the paper vs a single trial estimator:
print("NSW J=3 working-independence SE:", round(np.sqrt(1 / (1 / v[:3]).sum()), 1), "| DIM-RCT SE:", round(np.sqrt(v[0]), 1))
