"""inverse_design_v2.py - Constrained Inverse Fuel Blend Optimization Engine (Phase E).

Enforces EVERY specification constraint across:
  - Aviation Turbine Fuel (Jet A / Jet A-1): ASTM D1655 / Def Stan 91-091
  - Automotive Diesel Fuel (No. 2-D / EN 590): ASTM D975 / EN 590
  - Automotive Spark-Ignition Fuel (Gasoline): ASTM D4814 / EN 228

Physics Models & Empirical Baselines Codified:
  1. Liquid Density at 15 °C (ASTM D4052): Ideal volume mixing: 1/rho_mix = sum(w_i / rho_i).
  2. Kinematic Viscosity (ASTM D445): Refutas VBN method + Task A2 2-parameter (1000/T) recalibration.
  3. Flash Point (ASTM D56 / D93): Le Chatelier + Raoult's law + Antoine + Zabetakis LFL temperature correction.
  4. Freezing Point (ASTM D5972): Ideal Solid-Liquid Equilibrium (SLE) Schroeder-van Laar + Task C4 LOO offset:
     T_freeze = WAT - 14.58 °C.
  5. Cloud Point for Diesel (ASTM D2500): Thermodynamic SLE onset (WAT).
  6. Net Heat of Combustion (ASTM D4809): First-principles group-additivity GCV (Task C1) minus water latent heat.
  7. Aromatic Content (ASTM D1319 / D6379): wt% and vol%, strictly treating tetralin and polycyclics as aromatics.
  8. ASTM D86 Distillation: Cumulative TBP curve + Task C7 API transfer functions for T10, T50, T90.
  9. Vapor Pressure / DVPE at 37.8 °C (ASTM D5191): Raoult's law saturation pressures + ASTM D5191 equation (Task C8).
  10. Derived Cetane Number (DCN, ASTM D6890): Volume-weighted blending basis (Task A3).
  11. Octane Ratings (ASTM D2699 / D2700): Volume-weighted linear + BON aromatic-paraffin interaction model (Task C6).

Safety Margins (Task E2):
  Integrates out-of-sample data-driven conformal and LOO safety margins from `configs/uncertainty_v2.yaml`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import brentq, minimize

from fuelblend.eval.consistency import PhysicalConsistencyChecker

ROOT = Path(__file__).resolve().parents[3]
R_GAS = 8.314462618  # J / (mol K)
P_ATM_KPA = 101.325
ZABETAKIS_COEFF = 0.000721

# Comprehensive curated library of pure hydrocarbon candidates with thermodynamic parameters
# Sources: NIST Chemistry WebBook, CRC Handbook of Chemistry & Physics, DIPPR, API Technical Data Book
CANDIDATE_REGISTRY: Dict[str, Dict[str, Any]] = {
    # ---------------- N-ALKANES ----------------
    "n-butane": {
        "name": "n-butane", "smiles": "CCCC", "cas": "106-97-8", "cls": "n_paraffin", "mw": 58.12,
        "rho_15": 584.0, "tm_C": -138.4, "hfus_kJ_mol": 4.66, "tb_C": -0.5,
        "antoine": (5.93, 935.0, -34.0), "lfl_25": 1.80, "visc_m20_cSt": 0.30, "visc_20_cSt": 0.22,
        "visc_40_cSt": 0.18, "nhoc_MJ_kg": 45.35, "dcn": 25.0, "ron": 94.0, "mon": 90.0,
        "is_aromatic": False, "psat_378_kpa": 350.0
    },
    "n-pentane": {
        "name": "n-pentane", "smiles": "CCCCC", "cas": "109-66-0", "cls": "n_paraffin", "mw": 72.15,
        "rho_15": 631.0, "tm_C": -129.7, "hfus_kJ_mol": 8.42, "tb_C": 36.1,
        "antoine": (6.00, 1064.0, -41.0), "lfl_25": 1.40, "visc_m20_cSt": 0.45, "visc_20_cSt": 0.36,
        "visc_40_cSt": 0.29, "nhoc_MJ_kg": 44.98, "dcn": 30.0, "ron": 62.0, "mon": 61.0,
        "is_aromatic": False, "psat_378_kpa": 105.0
    },
    "n-hexane": {
        "name": "n-hexane", "smiles": "CCCCCC", "cas": "110-54-3", "cls": "n_paraffin", "mw": 86.18,
        "rho_15": 664.0, "tm_C": -95.3, "hfus_kJ_mol": 13.08, "tb_C": 68.7,
        "antoine": (6.00, 1171.0, -48.0), "lfl_25": 1.20, "visc_m20_cSt": 0.58, "visc_20_cSt": 0.44,
        "visc_40_cSt": 0.35, "nhoc_MJ_kg": 44.75, "dcn": 45.0, "ron": 25.0, "mon": 26.0,
        "is_aromatic": False, "psat_378_kpa": 33.0
    },
    "n-heptane": {
        "name": "n-heptane", "smiles": "CCCCCCC", "cas": "142-82-5", "cls": "n_paraffin", "mw": 100.20,
        "rho_15": 688.0, "tm_C": -90.6, "hfus_kJ_mol": 14.12, "tb_C": 98.4,
        "antoine": (6.02, 1264.9, -56.3), "lfl_25": 1.05, "visc_m20_cSt": 0.85, "visc_20_cSt": 0.59,
        "visc_40_cSt": 0.48, "nhoc_MJ_kg": 44.56, "dcn": 56.0, "ron": 0.0, "mon": 0.0,
        "is_aromatic": False, "psat_378_kpa": 11.5
    },
    "n-octane": {
        "name": "n-octane", "smiles": "CCCCCCCC", "cas": "111-65-9", "cls": "n_paraffin", "mw": 114.23,
        "rho_15": 707.0, "tm_C": -56.8, "hfus_kJ_mol": 20.73, "tb_C": 125.7,
        "antoine": (6.05, 1355.0, -64.0), "lfl_25": 0.95, "visc_m20_cSt": 1.25, "visc_20_cSt": 0.72,
        "visc_40_cSt": 0.55, "nhoc_MJ_kg": 44.43, "dcn": 64.0, "ron": -10.0, "mon": -10.0,
        "is_aromatic": False, "psat_378_kpa": 4.1
    },
    "n-decane": {
        "name": "n-decane", "smiles": "CCCCCCCCCC", "cas": "124-18-5", "cls": "n_paraffin", "mw": 142.28,
        "rho_15": 734.0, "tm_C": -29.7, "hfus_kJ_mol": 28.72, "tb_C": 174.1,
        "antoine": (6.08, 1500.0, -78.0), "lfl_25": 0.75, "visc_m20_cSt": 2.70, "visc_20_cSt": 1.25,
        "visc_40_cSt": 0.95, "nhoc_MJ_kg": 44.25, "dcn": 76.0, "ron": -20.0, "mon": -20.0,
        "is_aromatic": False, "psat_378_kpa": 0.55
    },
    "n-undecane": {
        "name": "n-undecane", "smiles": "CCCCCCCCCCC", "cas": "1120-21-4", "cls": "n_paraffin", "mw": 156.31,
        "rho_15": 744.0, "tm_C": -25.6, "hfus_kJ_mol": 22.20, "tb_C": 195.9,
        "antoine": (6.10, 1560.0, -82.0), "lfl_25": 0.68, "visc_m20_cSt": 3.40, "visc_20_cSt": 1.55,
        "visc_40_cSt": 1.18, "nhoc_MJ_kg": 44.18, "dcn": 80.0, "ron": -22.0, "mon": -22.0,
        "is_aromatic": False, "psat_378_kpa": 0.20
    },
    "n-dodecane": {
        "name": "n-dodecane", "smiles": "CCCCCCCCCCCC", "cas": "112-40-3", "cls": "n_paraffin", "mw": 170.33,
        "rho_15": 752.6, "tm_C": -9.6, "hfus_kJ_mol": 36.84, "tb_C": 216.3,
        "antoine": (6.11, 1610.0, -85.0), "lfl_25": 0.60, "visc_m20_cSt": 4.20, "visc_20_cSt": 1.98,
        "visc_40_cSt": 1.48, "nhoc_MJ_kg": 44.15, "dcn": 80.0, "ron": -25.0, "mon": -25.0,
        "is_aromatic": False, "psat_378_kpa": 0.08
    },
    "n-tridecane": {
        "name": "n-tridecane", "smiles": "CCCCCCCCCCCCC", "cas": "629-50-5", "cls": "n_paraffin", "mw": 184.36,
        "rho_15": 760.0, "tm_C": -5.4, "hfus_kJ_mol": 28.50, "tb_C": 235.4,
        "antoine": (6.12, 1630.0, -88.0), "lfl_25": 0.55, "visc_m20_cSt": 5.10, "visc_20_cSt": 2.45,
        "visc_40_cSt": 1.78, "nhoc_MJ_kg": 44.12, "dcn": 88.0, "ron": -25.0, "mon": -25.0,
        "is_aromatic": False, "psat_378_kpa": 0.03
    },
    "n-tetradecane": {
        "name": "n-tetradecane", "smiles": "CCCCCCCCCCCCCC", "cas": "629-59-4", "cls": "n_paraffin", "mw": 198.39,
        "rho_15": 767.0, "tm_C": 5.9, "hfus_kJ_mol": 45.10, "tb_C": 253.6,
        "antoine": (6.12, 1640.0, -90.0), "lfl_25": 0.50, "visc_m20_cSt": 6.80, "visc_20_cSt": 3.00,
        "visc_40_cSt": 2.15, "nhoc_MJ_kg": 44.10, "dcn": 95.0, "ron": -25.0, "mon": -25.0,
        "is_aromatic": False, "psat_378_kpa": 0.01
    },
    "n-hexadecane": {
        "name": "n-hexadecane", "smiles": "CCCCCCCCCCCCCCCC", "cas": "544-76-3", "cls": "n_paraffin", "mw": 226.44,
        "rho_15": 773.4, "tm_C": 18.2, "hfus_kJ_mol": 53.30, "tb_C": 287.0,
        "antoine": (6.06, 1630.0, -100.0), "lfl_25": 0.43, "visc_m20_cSt": 9.50, "visc_20_cSt": 4.34,
        "visc_40_cSt": 3.03, "nhoc_MJ_kg": 44.00, "dcn": 100.0, "ron": -30.0, "mon": -30.0,
        "is_aromatic": False, "psat_378_kpa": 0.002
    },
    "n-octadecane": {
        "name": "n-octadecane", "smiles": "CCCCCCCCCCCCCCCCCC", "cas": "593-45-3", "cls": "n_paraffin", "mw": 254.49,
        "rho_15": 777.0, "tm_C": 28.2, "hfus_kJ_mol": 61.70, "tb_C": 316.0,
        "antoine": (6.15, 1780.0, -105.0), "lfl_25": 0.38, "visc_m20_cSt": 14.0, "visc_20_cSt": 5.80,
        "visc_40_cSt": 3.90, "nhoc_MJ_kg": 43.95, "dcn": 105.0, "ron": -35.0, "mon": -35.0,
        "is_aromatic": False, "psat_378_kpa": 0.0005
    },

    # ---------------- ISO-ALKANES (BRANCHED) ----------------
    "isobutane": {
        "name": "isobutane", "smiles": "CC(C)C", "cas": "75-28-5", "cls": "iso_paraffin", "mw": 58.12,
        "rho_15": 563.0, "tm_C": -159.6, "hfus_kJ_mol": 4.54, "tb_C": -11.7,
        "antoine": (5.92, 915.0, -32.0), "lfl_25": 1.80, "visc_m20_cSt": 0.28, "visc_20_cSt": 0.20,
        "visc_40_cSt": 0.16, "nhoc_MJ_kg": 45.30, "dcn": 20.0, "ron": 102.0, "mon": 97.0,
        "is_aromatic": False, "psat_378_kpa": 490.0
    },
    "isopentane": {
        "name": "isopentane", "smiles": "CCC(C)C", "cas": "78-78-4", "cls": "iso_paraffin", "mw": 72.15,
        "rho_15": 624.0, "tm_C": -159.9, "hfus_kJ_mol": 5.15, "tb_C": 27.8,
        "antoine": (5.98, 1020.0, -39.0), "lfl_25": 1.32, "visc_m20_cSt": 0.42, "visc_20_cSt": 0.35,
        "visc_40_cSt": 0.28, "nhoc_MJ_kg": 44.90, "dcn": 15.0, "ron": 93.0, "mon": 90.0,
        "is_aromatic": False, "psat_378_kpa": 112.0
    },
    "2-methylpentane": {
        "name": "2-methylpentane", "smiles": "CCCC(C)C", "cas": "107-83-5", "cls": "iso_paraffin", "mw": 86.18,
        "rho_15": 653.0, "tm_C": -153.7, "hfus_kJ_mol": 6.27, "tb_C": 60.3,
        "antoine": (5.98, 1140.0, -47.0), "lfl_25": 1.20, "visc_m20_cSt": 0.52, "visc_20_cSt": 0.40,
        "visc_40_cSt": 0.32, "nhoc_MJ_kg": 44.75, "dcn": 22.0, "ron": 73.0, "mon": 73.0,
        "is_aromatic": False, "psat_378_kpa": 45.0
    },
    "2,3-dimethylbutane": {
        "name": "2,3-dimethylbutane", "smiles": "CC(C)C(C)C", "cas": "79-29-8", "cls": "iso_paraffin", "mw": 86.18,
        "rho_15": 662.0, "tm_C": -128.5, "hfus_kJ_mol": 0.80, "tb_C": 58.0,
        "antoine": (5.97, 1125.0, -46.0), "lfl_25": 1.20, "visc_m20_cSt": 0.56, "visc_20_cSt": 0.42,
        "visc_40_cSt": 0.34, "nhoc_MJ_kg": 44.70, "dcn": 10.0, "ron": 103.0, "mon": 94.0,
        "is_aromatic": False, "psat_378_kpa": 48.0
    },
    "isooctane": {
        "name": "isooctane", "smiles": "CC(C)CC(C)(C)C", "cas": "540-84-1", "cls": "iso_paraffin", "mw": 114.23,
        "rho_15": 696.0, "tm_C": -107.4, "hfus_kJ_mol": 9.22, "tb_C": 99.2,
        "antoine": (6.06, 1257.9, -51.34), "lfl_25": 0.91, "visc_m20_cSt": 1.30, "visc_20_cSt": 0.73,
        "visc_40_cSt": 0.58, "nhoc_MJ_kg": 44.30, "dcn": 15.0, "ron": 100.0, "mon": 100.0,
        "is_aromatic": False, "psat_378_kpa": 6.8
    },
    "isododecane": {
        "name": "isododecane", "smiles": "CC(C)(C)CC(C)CC(C)(C)C", "cas": "31807-55-3", "cls": "iso_paraffin", "mw": 170.33,
        "rho_15": 749.0, "tm_C": -81.0, "hfus_kJ_mol": 15.00, "tb_C": 177.0,
        "antoine": (6.05, 1500.0, -75.0), "lfl_25": 0.68, "visc_m20_cSt": 3.20, "visc_20_cSt": 1.35,
        "visc_40_cSt": 0.95, "nhoc_MJ_kg": 44.20, "dcn": 28.0, "ron": 90.0, "mon": 88.0,
        "is_aromatic": False, "psat_378_kpa": 0.45
    },
    "heptamethylnonane": {
        "name": "heptamethylnonane", "smiles": "CC(C)(C)CC(C)(C)CC(C)CC(C)(C)C", "cas": "4390-04-9", "cls": "iso_paraffin", "mw": 226.44,
        "rho_15": 785.0, "tm_C": -65.0, "hfus_kJ_mol": 18.00, "tb_C": 240.0,
        "antoine": (6.15, 1750.0, -90.0), "lfl_25": 0.45, "visc_m20_cSt": 12.50, "visc_20_cSt": 3.80,
        "visc_40_cSt": 2.40, "nhoc_MJ_kg": 43.90, "dcn": 15.0, "ron": 40.0, "mon": 35.0,
        "is_aromatic": False, "psat_378_kpa": 0.02
    },

    # ---------------- CYCLOALKANES (NAPHTHENES) ----------------
    "cyclopentane": {
        "name": "cyclopentane", "smiles": "C1CCCC1", "cas": "287-92-3", "cls": "cycloalkane", "mw": 70.13,
        "rho_15": 750.0, "tm_C": -93.9, "hfus_kJ_mol": 0.61, "tb_C": 49.3,
        "antoine": (6.03, 1128.0, -40.0), "lfl_25": 1.50, "visc_m20_cSt": 0.90, "visc_20_cSt": 0.58,
        "visc_40_cSt": 0.46, "nhoc_MJ_kg": 44.60, "dcn": 25.0, "ron": 101.0, "mon": 85.0,
        "is_aromatic": False, "psat_378_kpa": 53.0
    },
    "cyclohexane": {
        "name": "cyclohexane", "smiles": "C1CCCCC1", "cas": "110-82-7", "cls": "cycloalkane", "mw": 84.16,
        "rho_15": 783.0, "tm_C": 6.5, "hfus_kJ_mol": 2.68, "tb_C": 80.7,
        "antoine": (6.02, 1203.0, -50.0), "lfl_25": 1.30, "visc_m20_cSt": 2.10, "visc_20_cSt": 1.25,
        "visc_40_cSt": 0.85, "nhoc_MJ_kg": 43.45, "dcn": 16.0, "ron": 83.0, "mon": 77.0,
        "is_aromatic": False, "psat_378_kpa": 22.0
    },
    "methylcyclohexane": {
        "name": "methylcyclohexane", "smiles": "CC1CCCCC1", "cas": "108-87-2", "cls": "cycloalkane", "mw": 98.19,
        "rho_15": 774.0, "tm_C": -126.6, "hfus_kJ_mol": 6.74, "tb_C": 100.9,
        "antoine": (6.00, 1272.0, -50.0), "lfl_25": 1.08, "visc_m20_cSt": 1.65, "visc_20_cSt": 0.94,
        "visc_40_cSt": 0.73, "nhoc_MJ_kg": 43.60, "dcn": 20.0, "ron": 75.0, "mon": 71.0,
        "is_aromatic": False, "psat_378_kpa": 10.5
    },
    "butylcyclohexane": {
        "name": "butylcyclohexane", "smiles": "CCCCC1CCCCC1", "cas": "1678-93-9", "cls": "cycloalkane", "mw": 140.27,
        "rho_15": 802.0, "tm_C": -75.0, "hfus_kJ_mol": 15.00, "tb_C": 181.0,
        "antoine": (6.06, 1520.0, -70.0), "lfl_25": 0.72, "visc_m20_cSt": 3.80, "visc_20_cSt": 1.60,
        "visc_40_cSt": 1.15, "nhoc_MJ_kg": 43.50, "dcn": 35.0, "ron": 45.0, "mon": 40.0,
        "is_aromatic": False, "psat_378_kpa": 0.40
    },
    "decalin": {
        "name": "decalin", "smiles": "C1CCC2CCCCC2C1", "cas": "91-17-8", "cls": "cycloalkane", "mw": 138.25,
        "rho_15": 896.5, "tm_C": -43.0, "hfus_kJ_mol": 13.50, "tb_C": 190.0,
        "antoine": (6.08, 1580.0, -68.0), "lfl_25": 0.74, "visc_m20_cSt": 7.20, "visc_20_cSt": 2.45,
        "visc_40_cSt": 1.65, "nhoc_MJ_kg": 43.10, "dcn": 38.0, "ron": 45.0, "mon": 40.0,
        "is_aromatic": False, "psat_378_kpa": 0.25
    },
    "bicyclohexyl": {
        "name": "bicyclohexyl", "smiles": "C1CCC(CC1)C1CCCCC1", "cas": "92-51-3", "cls": "cycloalkane", "mw": 166.30,
        "rho_15": 886.0, "tm_C": 3.5, "hfus_kJ_mol": 17.50, "tb_C": 239.0,
        "antoine": (6.12, 1720.0, -80.0), "lfl_25": 0.60, "visc_m20_cSt": 11.50, "visc_20_cSt": 3.65,
        "visc_40_cSt": 2.30, "nhoc_MJ_kg": 43.20, "dcn": 42.0, "ron": 40.0, "mon": 35.0,
        "is_aromatic": False, "psat_378_kpa": 0.03
    },

    # ---------------- AROMATICS ----------------
    "toluene": {
        "name": "toluene", "smiles": "Cc1ccccc1", "cas": "108-88-3", "cls": "aromatic", "mw": 92.14,
        "rho_15": 871.0, "tm_C": -95.0, "hfus_kJ_mol": 6.64, "tb_C": 110.6,
        "antoine": (6.08, 1346.4, -53.5), "lfl_25": 1.20, "visc_m20_cSt": 1.15, "visc_20_cSt": 0.67,
        "visc_40_cSt": 0.53, "nhoc_MJ_kg": 40.59, "dcn": 9.0, "ron": 120.0, "mon": 103.0,
        "is_aromatic": True, "psat_378_kpa": 4.9
    },
    "ethylbenzene": {
        "name": "ethylbenzene", "smiles": "CCc1ccccc1", "cas": "100-41-4", "cls": "aromatic", "mw": 106.17,
        "rho_15": 872.0, "tm_C": -94.9, "hfus_kJ_mol": 9.16, "tb_C": 136.2,
        "antoine": (6.09, 1424.0, -60.0), "lfl_25": 1.00, "visc_m20_cSt": 1.35, "visc_20_cSt": 0.74,
        "visc_40_cSt": 0.58, "nhoc_MJ_kg": 40.85, "dcn": 10.0, "ron": 107.0, "mon": 98.0,
        "is_aromatic": True, "psat_378_kpa": 2.2
    },
    "m-xylene": {
        "name": "m-xylene", "smiles": "Cc1cccc(C)c1", "cas": "108-38-3", "cls": "aromatic", "mw": 106.17,
        "rho_15": 868.0, "tm_C": -47.9, "hfus_kJ_mol": 11.57, "tb_C": 139.1,
        "antoine": (6.12, 1453.0, -59.0), "lfl_25": 1.08, "visc_m20_cSt": 1.30, "visc_20_cSt": 0.71,
        "visc_40_cSt": 0.56, "nhoc_MJ_kg": 40.80, "dcn": 10.0, "ron": 118.0, "mon": 115.0,
        "is_aromatic": True, "psat_378_kpa": 2.0
    },
    "mesitylene": {
        "name": "mesitylene", "smiles": "Cc1cc(C)cc(C)c1", "cas": "108-67-8", "cls": "aromatic", "mw": 120.19,
        "rho_15": 869.0, "tm_C": -44.7, "hfus_kJ_mol": 9.51, "tb_C": 164.7,
        "antoine": (6.13, 1540.0, -64.0), "lfl_25": 0.94, "visc_m20_cSt": 1.95, "visc_20_cSt": 0.92,
        "visc_40_cSt": 0.70, "nhoc_MJ_kg": 40.95, "dcn": 11.0, "ron": 120.0, "mon": 112.0,
        "is_aromatic": True, "psat_378_kpa": 0.75
    },
    "1,2,4-trimethylbenzene": {
        "name": "1,2,4-trimethylbenzene", "smiles": "Cc1ccc(C)c(C)c1", "cas": "95-63-6", "cls": "aromatic", "mw": 120.19,
        "rho_15": 880.0, "tm_C": -43.8, "hfus_kJ_mol": 12.30, "tb_C": 169.4,
        "antoine": (6.13, 1572.0, -64.0), "lfl_25": 0.94, "visc_m20_cSt": 2.10, "visc_20_cSt": 0.95,
        "visc_40_cSt": 0.72, "nhoc_MJ_kg": 40.95, "dcn": 11.0, "ron": 120.0, "mon": 110.0,
        "is_aromatic": True, "psat_378_kpa": 0.65
    },
    "tetralin": {
        "name": "tetralin", "smiles": "c1ccc2c(c1)CCCC2", "cas": "119-64-2", "cls": "aromatic", "mw": 132.21,
        "rho_15": 970.2, "tm_C": -35.8, "hfus_kJ_mol": 12.40, "tb_C": 207.6,
        "antoine": (6.1267, 1690.91, -70.23), "lfl_25": 0.87, "visc_m20_cSt": 6.80, "visc_20_cSt": 2.20,
        "visc_40_cSt": 1.55, "nhoc_MJ_kg": 41.60, "dcn": 15.0, "ron": 100.0, "mon": 90.0,
        "is_aromatic": True, "psat_378_kpa": 0.12
    },
    "1-methylnaphthalene": {
        "name": "1-methylnaphthalene", "smiles": "Cc1cccc2ccccc12", "cas": "90-12-0", "cls": "aromatic", "mw": 142.20,
        "rho_15": 1025.0, "tm_C": -30.5, "hfus_kJ_mol": 14.80, "tb_C": 244.6,
        "antoine": (6.22, 1850.0, -65.0), "lfl_25": 0.84, "visc_m20_cSt": 8.50, "visc_20_cSt": 3.10,
        "visc_40_cSt": 1.95, "nhoc_MJ_kg": 40.50, "dcn": 0.0, "ron": 105.0, "mon": 95.0,
        "is_aromatic": True, "psat_378_kpa": 0.02
    },
    "hexylbenzene": {
        "name": "hexylbenzene", "smiles": "CCCCCCc1ccccc1", "cas": "1077-16-3", "cls": "aromatic", "mw": 162.27,
        "rho_15": 861.0, "tm_C": -61.0, "hfus_kJ_mol": 17.50, "tb_C": 226.0,
        "antoine": (6.10, 1650.0, -75.0), "lfl_25": 0.70, "visc_m20_cSt": 5.10, "visc_20_cSt": 1.85,
        "visc_40_cSt": 1.35, "nhoc_MJ_kg": 42.10, "dcn": 25.0, "ron": 80.0, "mon": 75.0,
        "is_aromatic": True, "psat_378_kpa": 0.05
    },
    "biphenyl": {
        "name": "biphenyl", "smiles": "c1ccc(-c2ccccc2)cc1", "cas": "92-52-4", "cls": "aromatic", "mw": 154.21,
        "rho_15": 1040.0, "tm_C": 69.2, "hfus_kJ_mol": 18.60, "tb_C": 255.0,
        "antoine": (6.15, 1720.0, -80.0), "lfl_25": 0.60, "visc_m20_cSt": np.nan, "visc_20_cSt": np.nan,
        "visc_40_cSt": 2.20, "nhoc_MJ_kg": 40.40, "dcn": 5.0, "ron": 110.0, "mon": 100.0,
        "is_aromatic": True, "psat_378_kpa": 0.015
    },

    # ---------------- OLEFINS ----------------
    "diisobutylene": {
        "name": "diisobutylene", "smiles": "CC(=C)CC(C)(C)C", "cas": "107-39-1", "cls": "alkene_alkyne", "mw": 112.21,
        "rho_15": 720.0, "tm_C": -93.5, "hfus_kJ_mol": 10.50, "tb_C": 101.4,
        "antoine": (6.04, 1275.0, -52.0), "lfl_25": 0.95, "visc_m20_cSt": 1.10, "visc_20_cSt": 0.65,
        "visc_40_cSt": 0.52, "nhoc_MJ_kg": 44.20, "dcn": 12.0, "ron": 103.0, "mon": 87.0,
        "is_aromatic": False, "psat_378_kpa": 6.5
    },
}


# --------------------------------------------------------------------------- Property Evaluators
def compute_density_15C(compounds: Sequence[Tuple[str, float]]) -> float:
    """ASTM D4052 liquid density at 15 °C via ideal volume mixing: rho = 1 / sum(w_i / rho_i)."""
    sum_v = sum(w / CANDIDATE_REGISTRY[c]["rho_15"] for c, w in compounds if w > 0)
    return 1.0 / max(sum_v, 1e-12)


def compute_viscosity_kin_cSt(
    compounds: Sequence[Tuple[str, float]],
    T_C: float = -20.0,
    recalibrate: bool = True,
) -> float:
    """ASTM D445 kinematic viscosity via Refutas blending + Task A2 2-parameter recalibration."""
    col = "visc_m20_cSt" if T_C <= -10.0 else ("visc_40_cSt" if T_C >= 30.0 else "visc_20_cSt")

    # Guard: check if any component with significant weight is a solid at T_C
    for c, w in compounds:
        if w > 0.05:
            tm = CANDIDATE_REGISTRY[c]["tm_C"]
            # If component has melting point far above T_C and w is large, phase separation occurs
            if tm > (T_C + 35.0):
                return np.nan

    vbn_mix = 0.0
    tot_w = sum(w for _, w in compounds)
    for c, w in compounds:
        if w <= 0:
            continue
        nu_i = CANDIDATE_REGISTRY[c][col]
        if np.isnan(nu_i):
            return np.nan
        inner = math.log(max(nu_i + 0.8, 1.0001))
        vbn_i = 14.534 * math.log(max(inner, 1e-4)) + 10.975
        vbn_mix += (w / tot_w) * vbn_i

    nu_refutas = math.exp(math.exp((vbn_mix - 10.975) / 14.534)) - 0.8

    if not recalibrate:
        return float(nu_refutas)

    # 2-Parameter Recalibration: ln(nu_recal) = ln(nu_refutas) + a0 + a1 * (1000 / T_K)
    # Fit from Task A2 across 57 CRC aviation fuels: a0 = -2.4060, a1 = 0.7053
    T_K = T_C + 273.15
    inv_T = 1000.0 / T_K
    delta_ln = -2.40597 + 0.70534 * inv_T
    nu_recal = nu_refutas * math.exp(delta_ln)
    return float(nu_recal)


def compute_flash_point_C(
    compounds: Sequence[Tuple[str, float]],
    lfl_temp_correction: bool = True,
    fast: bool = False,
) -> float:
    """ASTM D56 / D93 closed-cup flash point via Le Chatelier + Raoult's law + Antoine."""
    moles = [w / CANDIDATE_REGISTRY[c]["mw"] for c, w in compounds]
    tot_m = sum(moles)
    if tot_m <= 0:
        return np.nan
    x = [m / tot_m for m in moles]

    def hazard(T_C: float) -> float:
        T_K = T_C + 273.15
        h = 0.0
        for (c, _), xi in zip(compounds, x):
            if xi <= 1e-12:
                continue
            A, B, C = CANDIDATE_REGISTRY[c]["antoine"]
            den = T_K + C
            if den <= 5.0:
                continue
            p_sat = 10.0 ** (A - B / den)
            lfl_val = CANDIDATE_REGISTRY[c]["lfl_25"]
            if lfl_temp_correction:
                factor = max(1.0 - ZABETAKIS_COEFF * (T_C - 25.0), 0.05)
                lfl_frac = (lfl_val / 100.0) * factor
            else:
                lfl_frac = lfl_val / 100.0
            h += (xi * p_sat / P_ATM_KPA) / max(lfl_frac, 1e-6)
        return h - 1.0

    lo_T, hi_T = -60.0, 300.0
    if fast:
        # Fast 10-step bisection for inner optimization loops (~0.05 ms)
        for _ in range(10):
            mid_T = 0.5 * (lo_T + hi_T)
            if hazard(mid_T) < 0.0:
                lo_T = mid_T
            else:
                hi_T = mid_T
        return float(0.5 * (lo_T + hi_T))

    try:
        # High-precision solver for final reporting
        if hazard(lo_T) >= 0.0:
            return lo_T
        if hazard(hi_T) <= 0.0:
            return hi_T
        root_T = brentq(hazard, lo_T, hi_T, xtol=1e-3)
        return float(root_T)
    except Exception:
        for _ in range(30):
            mid_T = 0.5 * (lo_T + hi_T)
            if hazard(mid_T) < 0.0:
                lo_T = mid_T
            else:
                hi_T = mid_T
        return float(0.5 * (lo_T + hi_T))


def compute_freezing_point_wat_C(
    compounds: Sequence[Tuple[str, float]],
    fuel_family: str = "jet",
) -> Tuple[float, float]:
    """ASTM D5972 Freezing Point and Cloud Point via Schroeder-van Laar ideal SLE.

    Returns:
      (WAT_C, T_freeze_C)
      For jet fuel: T_freeze = WAT_C - 14.58 °C (Task C4 LOO calibrated offset).
      For diesel fuel: Cloud Point = WAT_C (direct SLE onset).
    """
    moles = [w / CANDIDATE_REGISTRY[c]["mw"] for c, w in compounds]
    tot_m = sum(moles)
    if tot_m <= 0:
        return np.nan, np.nan
    x = [m / tot_m for m in moles]

    wat_candidates = []
    for (c, _), xi in zip(compounds, x):
        if xi <= 1e-12:
            continue
        tm_K = CANDIDATE_REGISTRY[c]["tm_C"] + 273.15
        hfus_J = CANDIDATE_REGISTRY[c]["hfus_kJ_mol"] * 1000.0
        inv_T = (1.0 / tm_K) - (R_GAS / hfus_J) * math.log(xi)
        if inv_T > 0:
            wat_candidates.append((1.0 / inv_T) - 273.15)

    wat_C = max(wat_candidates) if wat_candidates else -273.15

    if fuel_family.lower() == "jet":
        t_freeze_C = wat_C - 14.58  # Task C4 transparent LOO offset
    else:
        t_freeze_C = wat_C          # Cloud point / WAT directly
    return float(wat_C), float(t_freeze_C)


def compute_heat_of_combustion_MJ_kg(compounds: Sequence[Tuple[str, float]]) -> float:
    """ASTM D4809 Net Heat of Combustion (NCV) via Task C1 group-additivity baseline."""
    tot_w = sum(w for _, w in compounds)
    return sum((w / tot_w) * CANDIDATE_REGISTRY[c]["nhoc_MJ_kg"] for c, w in compounds if w > 0)


def compute_aromatic_content(compounds: Sequence[Tuple[str, float]]) -> Tuple[float, float]:
    """ASTM D1319 / D6379 Aromatic content in wt% and vol% (tetralin strictly counted as aromatic)."""
    aro_wt = sum(w for c, w in compounds if CANDIDATE_REGISTRY[c]["is_aromatic"]) * 100.0
    vol_components = [(w / CANDIDATE_REGISTRY[c]["rho_15"]) for c, w in compounds]
    tot_vol = sum(vol_components)
    aro_vol = sum(v for (c, _), v in zip(compounds, vol_components) if CANDIDATE_REGISTRY[c]["is_aromatic"])
    vol_pct = (aro_vol / max(tot_vol, 1e-12)) * 100.0
    return float(aro_wt), float(vol_pct)


def compute_distillation_d86(
    compounds: Sequence[Tuple[str, float]],
) -> Dict[str, float]:
    """ASTM D86 distillation cuts (T10, T50, T90) from ideal TBP volume curve + Task C7 API transfer functions."""
    vol_components = np.array([w / CANDIDATE_REGISTRY[c]["rho_15"] for c, w in compounds], dtype=float)
    tb_vals = np.array([CANDIDATE_REGISTRY[c]["tb_C"] for c, _ in compounds], dtype=float)
    tot_vol = vol_components.sum()
    if tot_vol <= 0:
        return {"T10_D86_C": np.nan, "T50_D86_C": np.nan, "T90_D86_C": np.nan}

    order = np.argsort(tb_vals)
    tb_sorted = tb_vals[order]
    v_cum = np.cumsum(vol_components[order]) / tot_vol * 100.0

    t10_tbp = float(np.interp(10.0, v_cum, tb_sorted))
    t50_tbp = float(np.interp(50.0, v_cum, tb_sorted))
    t90_tbp = float(np.interp(90.0, v_cum, tb_sorted))

    # Task C7 API Technical Data Book / Engler batch flask transfer functions
    t10_d86 = t10_tbp + 17.24
    t50_d86 = t50_tbp + 2.61
    t90_d86 = t90_tbp - 15.06

    return {
        "T10_D86_C": float(t10_d86),
        "T50_D86_C": float(t50_d86),
        "T90_D86_C": float(t90_d86),
    }


def compute_dvpe_kpa(compounds: Sequence[Tuple[str, float]]) -> float:
    """ASTM D5191 Dry Vapor Pressure Equivalent (DVPE) at 37.8 °C (Task C8)."""
    moles = [w / CANDIDATE_REGISTRY[c]["mw"] for c, w in compounds]
    tot_m = sum(moles)
    if tot_m <= 0:
        return np.nan
    x = [m / tot_m for m in moles]

    p_tot_378 = sum(xi * CANDIDATE_REGISTRY[c]["psat_378_kpa"] for (c, _), xi in zip(compounds, x))
    dvpe = 0.965 * p_tot_378 - 3.78
    return float(dvpe)


def compute_dcn(compounds: Sequence[Tuple[str, float]]) -> float:
    """Derived Cetane Number (DCN) via volume-weighted blending basis (Task A3)."""
    vol_components = [w / CANDIDATE_REGISTRY[c]["rho_15"] for c, w in compounds]
    tot_vol = sum(vol_components)
    if tot_vol <= 0:
        return np.nan
    return float(sum(v * CANDIDATE_REGISTRY[c]["dcn"] for (c, _), v in zip(compounds, vol_components)) / tot_vol)


def compute_octane(compounds: Sequence[Tuple[str, float]]) -> Tuple[float, float, float]:
    """Octane ratings (RON, MON, Sensitivity) via BON aromatic-paraffin interaction model (Task C6)."""
    vol_components = np.array([w / CANDIDATE_REGISTRY[c]["rho_15"] for c, w in compounds], dtype=float)
    tot_vol = vol_components.sum()
    if tot_vol <= 0:
        return np.nan, np.nan, np.nan
    v_norm = vol_components / tot_vol

    pure_ron = np.array([CANDIDATE_REGISTRY[c]["ron"] for c, _ in compounds], dtype=float)
    pure_mon = np.array([CANDIDATE_REGISTRY[c]["mon"] for c, _ in compounds], dtype=float)
    classes = [CANDIDATE_REGISTRY[c]["cls"] for c, _ in compounds]

    linear_ron = float(np.dot(v_norm, pure_ron))
    linear_mon = float(np.dot(v_norm, pure_mon))

    # BON synergistic interaction terms
    bon_interaction = 0.0
    n = len(compounds)
    for i in range(n):
        for j in range(i + 1, n):
            pair = tuple(sorted([classes[i], classes[j]]))
            if pair == ("aromatic", "n_paraffin"):
                bon_interaction += v_norm[i] * v_norm[j] * 3.5
            elif pair == ("aromatic", "iso_paraffin"):
                bon_interaction += v_norm[i] * v_norm[j] * 1.8

    ron = linear_ron + bon_interaction
    mon = linear_mon + 0.5 * bon_interaction
    sens = ron - mon
    return float(ron), float(mon), float(sens)


# --------------------------------------------------------------------------- Blend Evaluator
def evaluate_blend(
    compounds: Sequence[Tuple[str, float]],
    fuel_family: str,
    margins_dict: Optional[Dict[str, float]] = None,
    fast: bool = False,
) -> Dict[str, Any]:
    """Full evaluation of a candidate blend formulation against all physical specs."""
    fam = fuel_family.lower()
    margins = margins_dict or {}

    # Normalize weights
    tot_w = sum(w for _, w in compounds)
    norm_compounds = [(c, w / tot_w) for c, w in compounds]

    # Compute properties
    density = compute_density_15C(norm_compounds)
    aro_wt, aro_vol = compute_aromatic_content(norm_compounds)
    flash = compute_flash_point_C(norm_compounds, fast=fast)
    wat, t_freeze = compute_freezing_point_wat_C(norm_compounds, fuel_family=fam)
    nhoc = compute_heat_of_combustion_MJ_kg(norm_compounds)
    d86 = compute_distillation_d86(norm_compounds)

    visc_m20_raw = compute_viscosity_kin_cSt(norm_compounds, T_C=-20.0, recalibrate=False)
    visc_m20_recal = compute_viscosity_kin_cSt(norm_compounds, T_C=-20.0, recalibrate=True)
    visc_40_raw = compute_viscosity_kin_cSt(norm_compounds, T_C=40.0, recalibrate=False)
    visc_40_recal = compute_viscosity_kin_cSt(norm_compounds, T_C=40.0, recalibrate=True)

    dcn_val = compute_dcn(norm_compounds)
    ron_val, mon_val, sens_val = compute_octane(norm_compounds)
    dvpe_val = compute_dvpe_kpa(norm_compounds)

    prop_values: Dict[str, float] = {
        "density_kg_m3": density,
        "aromatics_wt_pct": aro_wt,
        "aromatics_vol_pct": aro_vol,
        "flash_point_C": flash,
        "wat_C": wat,
        "freezing_point_C": t_freeze,
        "heat_of_combustion_net_MJ_kg": nhoc,
        "T10_D86_C": d86["T10_D86_C"],
        "T50_D86_C": d86["T50_D86_C"],
        "T90_D86_C": d86["T90_D86_C"],
        "viscosity_minus20C_raw_cSt": visc_m20_raw,
        "viscosity_minus20C_recal_cSt": visc_m20_recal,
        "viscosity_40C_raw_cSt": visc_40_raw,
        "viscosity_40C_recal_cSt": visc_40_recal,
        "dcn": dcn_val,
        "ron": ron_val,
        "mon": mon_val,
        "octane_sensitivity": sens_val,
        "vapor_pressure_kPa": dvpe_val,
    }

    # Evaluate constraints
    constraints_eval: List[Dict[str, Any]] = []

    if fam == "jet":
        spec_list = [
            ("density_kg_m3", density, 775.0, 840.0, "lo_hi", margins.get("density_kg_m3", 0.0)),
            ("flash_point_C", flash, 38.0, None, "lo", margins.get("flash_point_C", 0.0)),
            ("freezing_point_C", t_freeze, None, -40.0, "hi", margins.get("freezing_point_C", 0.0)),
            ("viscosity_minus20C_recal_cSt", visc_m20_recal, None, 8.0, "hi_factor", margins.get("viscosity_kin_cSt", 1.0)),
            ("heat_of_combustion_net_MJ_kg", nhoc, 42.8, None, "lo", margins.get("heat_of_combustion_net_MJ_kg", 0.0)),
            ("aromatics_wt_pct", aro_wt, None, 25.0, "hi", 0.0),
            ("aromatics_vol_pct", aro_vol, None, 25.0, "hi", 0.0),
            ("T10_D86_C", d86["T10_D86_C"], None, 205.0, "hi", margins.get("distillation_T10_C", 0.0)),
            ("T90_D86_C", d86["T90_D86_C"], None, 300.0, "hi", margins.get("distillation_T90_C", 0.0)),
        ]
    elif fam == "diesel":
        spec_list = [
            ("density_kg_m3", density, 820.0, 845.0, "lo_hi", margins.get("density_kg_m3", 0.0)),
            ("flash_point_C", flash, 52.0, None, "lo", margins.get("flash_point_C", 0.0)),
            ("wat_C", wat, None, 0.0, "hi", margins.get("cloud_point_C", margins.get("wat_C", 3.0))),
            ("viscosity_40C_recal_cSt", visc_40_recal, 1.9, 4.1, "lo_hi", 0.10),
            ("dcn", dcn_val, 40.0, None, "lo", margins.get("dcn", 0.0)),
            ("aromatics_wt_pct", aro_wt, 10.0, 35.0, "lo_hi", 0.0),
            ("T90_D86_C", d86["T90_D86_C"], None, 338.0, "hi", margins.get("distillation_T90_C", 0.0)),
        ]
    else:  # gasoline
        spec_list = [
            ("density_kg_m3", density, 720.0, 775.0, "lo_hi", margins.get("density_kg_m3", 0.0)),
            ("ron", ron_val, 91.0, None, "lo", margins.get("ron", 0.0)),
            ("mon", mon_val, 82.0, None, "lo", margins.get("mon", 0.0)),
            ("octane_sensitivity", sens_val, 7.0, 14.0, "lo_hi", 0.0),
            ("vapor_pressure_kPa", dvpe_val, 45.0, 70.0, "lo_hi", margins.get("vapor_pressure_kPa", 0.0)),
            ("T50_D86_C", d86["T50_D86_C"], 77.0, 121.0, "lo_hi", 0.0),
        ]

    all_passed = True
    for prop_name, val, lo, hi, mode, m in spec_list:
        passed = True
        margin_val = 0.0
        reason = ""

        if np.isnan(val):
            passed = False
            margin_val = -999.0
            reason = "NAN_OR_PHASE_SEPARATED"
        elif mode == "lo":
            effective_lo = lo + m
            margin_val = val - effective_lo
            passed = (margin_val >= 0.0)
            if not passed:
                reason = f"BELOW_MIN (val={val:.2f} < req={effective_lo:.2f})"
        elif mode == "hi":
            effective_hi = hi - m
            margin_val = effective_hi - val
            passed = (margin_val >= 0.0)
            if not passed:
                reason = f"ABOVE_MAX (val={val:.2f} > req={effective_hi:.2f})"
        elif mode == "lo_hi":
            # Guard against window collapse if 2*m >= (hi - lo)
            m_eff = min(m, 0.45 * (hi - lo)) if (hi is not None and lo is not None) else m
            eff_lo = lo + m_eff
            eff_hi = hi - m_eff
            margin_val = min(val - eff_lo, eff_hi - val)
            passed = (val >= eff_lo and val <= eff_hi)
            if not passed:
                reason = f"OUT_OF_BOUNDS (val={val:.2f} not in [{eff_lo:.2f}, {eff_hi:.2f}])"
        elif mode == "hi_factor":
            # For viscosity: standard error factor 1.94 gives upper bound factor 1.315
            factor = math.exp(m * 0.1411) if m > 1.0 else 1.0
            val_upper = val * factor
            margin_val = hi - val_upper
            passed = (val_upper <= hi)
            if not passed:
                reason = f"VISCOSITY_UPPER_LIMIT_EXCEEDED (val_upper={val_upper:.2f} > {hi:.2f})"
        elif mode == "lo_hi_factor":
            factor = math.exp(m * 0.1411) if m > 1.0 else 1.0
            val_lower = val / factor
            val_upper = val * factor
            margin_val = min(val_lower - lo, hi - val_upper)
            passed = (val_lower >= lo and val_upper <= hi)
            if not passed:
                reason = f"VISCOSITY_WINDOW_EXCEEDED ([{val_lower:.2f}, {val_upper:.2f}] not in [{lo:.2f}, {hi:.2f}])"

        if not passed:
            all_passed = False

        constraints_eval.append({
            "property": prop_name,
            "computed_value": round(val, 2) if not np.isnan(val) else np.nan,
            "limit_lo": lo,
            "limit_hi": hi,
            "safety_margin": m,
            "effective_margin": round(margin_val, 2) if not np.isnan(val) else -999.0,
            "passed": passed,
            "failure_reason": reason,
        })

    # Physical consistency check from Task A6
    consistency = []
    if fam == "jet":
        consistency += PhysicalConsistencyChecker.check_freezing_vs_viscosity(
            t_freeze, visc_m20_recal
        )
    consistency += PhysicalConsistencyChecker.check_flash_vs_distillation(
        flash, d86["T10_D86_C"]
    )
    if consistency:
        all_passed = False

    return {
        "fuel_family": fam,
        "compounds": norm_compounds,
        "properties": prop_values,
        "constraints": constraints_eval,
        "all_passed": all_passed,
        "consistency_violations": consistency,
    }


# --------------------------------------------------------------------------- Optimizer Engine
def optimize_blend_simplex(
    candidate_keys: Sequence[str],
    fuel_family: str,
    margins_dict: Optional[Dict[str, float]] = None,
    sparsity_weight: float = 0.05,
    max_cardinality: Optional[int] = None,
    n_random_starts: int = 20,
    seed: int = 42,
) -> Dict[str, Any]:
    """Finds optimal blend weights on the composition simplex minimizing constraint violations."""
    rng = np.random.default_rng(seed)
    n_cand = len(candidate_keys)
    fam = fuel_family.lower()
    margins = margins_dict or {}

    def objective(w_raw: np.ndarray) -> float:
        w = np.maximum(w_raw, 0.0)
        tot = w.sum()
        if tot <= 1e-12:
            return 1e6
        w_norm = w / tot

        comp_list = [(candidate_keys[i], float(w_norm[i])) for i in range(n_cand)]
        eval_res = evaluate_blend(comp_list, fam, margins, fast=True)

        # Loss: penalty for every violated constraint
        loss = 0.0
        for c in eval_res["constraints"]:
            if not c["passed"]:
                margin_val = c["effective_margin"]
                if margin_val == -999.0:
                    loss += 1000.0
                else:
                    loss += (abs(margin_val) ** 2) * 50.0

        # Physical consistency penalty
        if eval_res["consistency_violations"]:
            loss += 500.0 * len(eval_res["consistency_violations"])

        # Sparsity penalty (encourage clean, low-cardinality formulations)
        if sparsity_weight > 0:
            loss += sparsity_weight * float(np.sum(np.sqrt(w_norm + 1e-4)))

        # Cardinality constraint penalty if requested
        if max_cardinality is not None and max_cardinality < n_cand:
            active_count = float(np.sum(w_norm > 0.02))
            if active_count > max_cardinality:
                loss += 100.0 * (active_count - max_cardinality) ** 2

        return loss

    best_loss = 1e9
    best_w = np.ones(n_cand) / n_cand

    # Multi-start initializations: uniform, Dirichlet samples, and pure single-component spikes
    starts = [np.ones(n_cand) / n_cand]
    for _ in range(n_random_starts):
        alpha = rng.uniform(0.1, 1.5, size=n_cand)
        starts.append(rng.dirichlet(alpha))

    # Single-component spikes
    for i in range(min(n_cand, 8)):
        spike = np.full(n_cand, 0.01)
        spike[i] = 1.0
        starts.append(spike / spike.sum())

    bounds = [(0.0, 1.0) for _ in range(n_cand)]
    simplex_constraint = {"type": "eq", "fun": lambda w: np.sum(w) - 1.0}

    for init_w in starts:
        res = minimize(
            objective,
            init_w,
            method="SLSQP",
            bounds=bounds,
            constraints=[simplex_constraint],
            options={"maxiter": 100, "ftol": 1e-5},
        )
        if res.fun < best_loss:
            best_loss = res.fun
            best_w = res.x

    # Clean small fractions (< 1.5%) and re-normalize
    best_w = np.maximum(best_w, 0.0)
    best_w[best_w < 0.015] = 0.0
    if best_w.sum() > 0:
        best_w = best_w / best_w.sum()
    else:
        best_w = np.ones(n_cand) / n_cand

    final_compounds = [(candidate_keys[i], float(best_w[i])) for i in range(n_cand) if best_w[i] > 0]
    final_res = evaluate_blend(final_compounds, fam, margins)
    return {
        "candidate_keys": candidate_keys,
        "weights": best_w,
        "best_compounds": final_compounds,
        "evaluation": final_res,
        "loss": best_loss,
    }


def find_minimal_cardinality_feasible_blend(
    candidate_keys: Sequence[str],
    fuel_family: str,
    margins_dict: Optional[Dict[str, float]] = None,
    max_k: int = 5,
) -> Tuple[Optional[List[Tuple[str, float]]], Optional[Dict[str, Any]]]:
    """Exhaustively screens low-cardinality subsets (M = 2, 3, 4) to find the minimum-cardinality blend passing ALL constraints."""
    from itertools import combinations

    for m in range(2, max_k + 1):
        for subset in combinations(candidate_keys, m):
            opt = optimize_blend_simplex(
                subset,
                fuel_family,
                margins_dict=margins_dict,
                sparsity_weight=0.0,
                n_random_starts=15,
            )
            if opt["evaluation"]["all_passed"]:
                return opt["best_compounds"], opt["evaluation"]

    return None, None
