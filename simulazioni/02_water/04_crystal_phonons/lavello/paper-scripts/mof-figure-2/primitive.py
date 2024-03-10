# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 15-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

import pathlib

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

d = pathlib.Path().absolute()
cifs = d.glob("*.cif")

for f in cifs:
   print(f"process {f}")
   stem = pathlib.Path(f).stem
   try:
     s = Structure.from_file(f, primitive=True)
     s.to(filename=f"prim/{stem}-primitive.cif")
     sga = SpacegroupAnalyzer(s,symprec=0.01)
     sc = sga.get_conventional_standard_structure()
     sc.to(filename=f"conv/{stem}-conventional.cif")
   except:
     continue


