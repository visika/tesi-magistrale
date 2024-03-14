# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 31-01-2024
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html
from __future__ import annotations

from typing import TYPE_CHECKING
import torch

architectures = ["mace", "mace_mp", "mace_off","m3gnet", "chgnet"]

if TYPE_CHECKING:
  from typing import LiteralString
  from ase.calculator.calculator import Calculator

"""
similar in spirit with matcalc and quacc approaches
- https://github.com/materialsvirtuallab/matcalc
- https://github.com/Quantum-Accelerators/quacc.git
"""

def choose_calculator( architecture: LiteralString[architectures] = 'mace', **kwargs ) -> Calculator:

  match architecture:
    case "mace":
      from mace import __version__
      from mace.calculators import MACECalculator

      if "default_dtype" not in kwargs:
        kwargs["default_dtype"] = "float64"
      if "device" not in kwargs:
        kwargs["device"] = "cuda"
      calculator = MACECalculator(**kwargs)

    case "mace_mp":
      from mace import __version__
      from mace.calculators import mace_mp

      if "default_dtype" not in kwargs:
        kwargs["default_dtype"] = "float64"
      if "device" not in kwargs:
        kwargs["device"] = "cuda"
      if "model" not in kwargs:
        kwargs["model"] = "small"
      calculator = mace_mp(**kwargs)

    case "mace_off":

      from mace import __version__
      from mace.calculators import mace_off

      if "default_dtype" not in kwargs:
        kwargs["default_dtype"] = "float64"
      if "device" not in kwargs:
        kwargs["device"] = "cuda"
      if "model" not in kwargs:
        kwargs["model"] = "small"
      calculator = mace_off(**kwargs)

    case "m3gnet":
      from matgl import __version__
      from matgl import load_model
      from matgl.ext.ase import M3GNetCalculator
      print(__version__)
      if "model" not in kwargs:
        model = load_model("M3GNet-MP-2021.2.8-DIRECT-PES")
      if "stress_weight" not in kwargs:
        kwargs.setdefault("stress_weight", 1.0 / 160.21766208)
      calculator = M3GNetCalculator(potential=model, **kwargs)

    case "chgnet":
      from chgnet import __version__
      from chgnet.model.dynamics import CHGNetCalculator
      if "use_device" not in kwargs:
        kwargs["use_device"] = "cuda"

      calculator = CHGNetCalculator(**kwargs)

    case _:
        raise ValueError(f"Unrecognized {architecture=}. suported architectures are {architectures}")

  calculator.parameters['version'] = __version__

  return calculator
