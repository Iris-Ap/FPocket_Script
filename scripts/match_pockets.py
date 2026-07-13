"""
match_pockets.py  --  Spatially match each fpocket-selected pocket to a DoGSite pocket,
and read the DoGSite enclosure/depth/drugScore of the match.

Matching = for each DoGSite pocket, count how many of the fpocket pocket's alpha
spheres fall within 4 A of that pocket's lining residues; the pocket with the most
overlap is the match. overlap == 0 => the two tools do NOT detect the same pocket
(a finding in itself, e.g. FOXM1).

Run:  python3 match_pockets.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

if __name__ == "__main__":
    print(f"{'protein':12} {'fpk':5} {'->DoGSite':10} {'overlap':10} "
          f"{'enclosure':10} {'depth':7} {'drugScore':9} {'quality'}")
    for name in C.ORDER:
        pk = C.PROTEINS[name]["fp_pocket"]
        pid, ov, nsp = C.match_dogsite_pocket(name, pk)
        d = C.dogsite_desc(name).get(pid, {})
        enc, dep, ds = d.get("enclosure", "-"), d.get("depth", "-"), d.get("drugScore", "-")
        frac = ov / nsp if nsp else 0
        quality = "none" if ov == 0 else ("strong" if frac > 0.3 else "weak")
        print(f"{name:12} P{pk:<4} {pid or '-':10} {str(ov)+'/'+str(nsp):10} "
              f"{enc:10} {dep:7} {ds:9} {quality}")
