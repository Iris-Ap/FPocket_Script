"""
pull_rubric_data.py  --  Extract the rubric metrics from the fpocket highest-druggability
pocket of every protein.

Rubric metrics (fpocket-native): Volume, Hydrophobicity, Polarity (% polar atoms).
Druggability Score is used ONLY to select the pocket, not as a rubric metric.
Enclosure is added separately via compute_enclosure.py.

Run:  python3 pull_rubric_data.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

FIELDS = [
    ("Volume (A^3)",        "Volume"),
    ("Hydrophobicity",      "Hydrophobicity score"),
    ("Polarity (%polar)",   "Proportion of polar atoms"),
]

if __name__ == "__main__":
    hdr = f"{'protein':12} {'type':8} {'pocket':7} " + " ".join(f"{n:18}" for n, _ in FIELDS)
    print(hdr)
    for name in C.ORDER:
        info = C.fpocket_info(name)
        pk = C.highest_druggability_pocket(name)   # = PROTEINS[name]['fp_pocket']
        d = info[pk]
        vals = " ".join(f"{float(d[f]):<18.2f}" for _, f in FIELDS)
        print(f"{name:12} {C.PROTEINS[name]['type']:8} P{pk:<6} {vals}")
