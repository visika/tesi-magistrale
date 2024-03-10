# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 31-01-2024
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

from mlip_calculators import choose_calculator

from ase import build

benzene = build.molecule('C6H6')

benzene.calc = choose_calculator(architecture="mace",device="cpu",model_paths="/home/drFaustroll/.cache/mace/5yyxdm76")
print(f"E_config= {benzene.get_potential_energy()} eV")
