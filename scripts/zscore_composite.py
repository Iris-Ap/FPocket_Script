"""
zscore_composite.py -- Build the weighted z-score druggability composite.

Baseline = the 3 druggable controls (ERα, HER2, PI3Kα): each metric's mean and
sample SD (n-1) come from them. Every protein is scored as z = (value - control mean) /
control SD, oriented so higher = more druggable (polarity flipped), then combined as a
weighted average of z-scores.

Raw values: Volume/Hydrophobicity/Polarity from fpocket (highest-druggability pocket);
Enclosure (burial) from grid_enclosure.py.

Run: python3 zscore_composite.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
import grid_enclosure as GE
from statistics import mean, stdev

METRICS = ["volume", "hydrophobicity", "polarity", "enclosure"]
WEIGHTS = {"volume": 40, "hydrophobicity": 90, "polarity": 60, "enclosure": 90}
FLIP    = {"polarity"}   # higher = LESS druggable -> flip z sign


def raw_values(name):
    d = C.fpocket_info(name)[C.PROTEINS[name]["fp_pocket"]]
    eo, _, _ = GE.grid_enclosure(name, C.PROTEINS[name]["fp_pocket"])
    return {
        "volume":        float(d["Volume"]) ** 0.5,   # sqrt(volume), following Halgren (2009)
        "hydrophobicity":float(d["Hydrophobicity score"]),
        "polarity":      float(d["Proportion of polar atoms"]),
        "enclosure":     1 - eo,          # burial
    }


raw = {n: raw_values(n) for n in C.ORDER}
controls = [n for n in C.ORDER if C.PROTEINS[n]["type"] == "control"]
stats = {m: (mean(raw[c][m] for c in controls), stdev(raw[c][m] for c in controls)) for m in METRICS}

def zscore(name, m):
    mu, sd = stats[m]
    z = (raw[name][m] - mu) / sd
    return -z if m in FLIP else z

Wsum = sum(WEIGHTS.values())
def composite(name):
    return sum(WEIGHTS[m] * zscore(name, m) for m in METRICS) / Wsum


if __name__ == "__main__":
    print("=== control baseline (mean, sample SD) per metric ===")
    for m in METRICS:
        mu, sd = stats[m]
        print(f"  {m:15} mean={mu:9.3f}  SD={sd:9.3f}")

    print("\n=== oriented z-scores (higher = more druggable) ===")
    hdr = f"{'protein':8}{'type':9}" + "".join(f"{m[:5]:>9}" for m in METRICS) + f"{'COMPOSITE':>12}"
    print(hdr)
    rows = []
    for n in C.ORDER:
        zs = [zscore(n, m) for m in METRICS]
        comp = composite(n)
        rows.append((n, comp))
        line = f"{n.split('_')[0]:8}{C.PROTEINS[n]['type']:9}" + "".join(f"{z:9.2f}" for z in zs) + f"{comp:12.3f}"
        print(line)

    print("\n=== ranked by composite (most -> least druggable) ===")
    for n, comp in sorted(rows, key=lambda r: -r[1]):
        print(f"  {n.split('_')[0]:8}{C.PROTEINS[n]['type']:9}{comp:8.3f}")
