# UCSF Chimera command script -- SRA_iris  (minimal, robust version)
# ERalpha protein (gray) + tamoxifen (red) + fpocket Pocket 1 (orange spheres).

open /Users/iris/SRA_iris/protein_clean.pdb
open /Users/iris/SRA_iris/drug_only.pdb
open /Users/iris/SRA_iris/pocket1.pdb

background solid white
~display
ribbon #0
color gray #0

display #1
represent stick #1
color red #1

display #2
represent sphere #2
color orange #2

focus
wait 30
copy file /Users/iris/SRA_iris/scene.png png supersample 3
