"""v4 Task 17: inverse-design engine that reads specs and safety margins from configuration ONLY.

Contract (fixes F5 of the v3 audit: v3 hard-coded the spec limits and typed the margins in the script):

* Specification limits come ONLY from ``configs/specs_<standard>.yaml`` (``families.<family>.constraints``).
  The ``uncertainty_margin*`` keys those files still carry (hand-typed in v3) are NEVER read.
* Safety margins come ONLY from ``configs/uncertainty_v4.yaml`` and ONLY through ``applied_margin``
  (never the raw conformal quantile in ``value``: a negative quantile is not extra safety).  An entry with no
  usable ``applied_margin`` (the generator writes ``value: inf`` and omits ``applied_margin``) means NO VALID
  MARGIN: the constraint cannot be certified and the engine reports ``uncertified`` with the reason.  No number
  is ever substituted.
* The flash-point margin of an optimizer-designed blend is the LARGER of the lab and the D93 ``applied_margin``
  (a model's ``margin_keys`` list is combined with ``max``), because the D93 margin floors at 0 (every CRC fuel is
  under-predicted) and 0 must not be used for a novel blend.
* A constraint status is one of
    ``pass``           finite margin from the yaml, margin-adjusted prediction satisfies the limit
    ``by_construction`` quantity fixed by the blend recipe itself (e.g. class mass fraction): no model error
    ``pass_jet_derived_margin_only``  meets the limit with the margin, but the margin was validated on jet fuels only and
                       this is a diesel/gasoline spec: NOT a certification, not counted as ok
    ``fail``           the (margin-adjusted, or nominal when no margin exists) prediction violates the limit
    ``uncertified``    cannot be certified: no valid margin, no model, unverified inputs, T/domain outside the
                       margin's validated scope, extrapolation, or a margin on the wrong side
  ``passed`` is True only when EVERY constraint is ``pass`` / ``by_construction`` and the design is not
  extrapolating.  A nominal violation is reported as ``fail`` even when no margin exists (a violated limit does
  not become acceptable because the margin is missing).
* ``extrapolating``: if more than ``EXTRAP_WT_THRESHOLD`` of the blend (by mass) belongs to classes absent from
  the validation data of a property model the design uses (or the pool contains no in-domain compound at all),
  the design is flagged and ``passed`` is False.

Models are injected (``models`` dict) so the engine itself holds no property physics; the physics wrappers live in
``fuelblend.design.models_v4`` and are built from the Task 4-15 library code.
"""
from __future__ import annotations

import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np
import yaml
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[3]
CONFIGS = ROOT / "configs"
UNCERTAINTY_YAML = CONFIGS / "uncertainty_v4.yaml"

# Non-spec engineering constants (the AST test allows exactly these numeric literals; see tests/test_v4_design.py).
EXTRAP_WT_THRESHOLD = 0.6      # wt fraction of absent-class mass above which a design is flagged extrapolating
BAD_SLACK = -100.0             # normalised slack assigned to a non-finite prediction inside the optimiser
MAX_ITER = 120

STATUS_PASS = "pass"
STATUS_CONSTRUCTION = "by_construction"
STATUS_FAIL = "fail"
STATUS_UNCERTIFIED = "uncertified"
STATUS_SCOPE = "pass_jet_derived_margin_only"   # meets the limit with the margin, but the margin was not validated for this fuel kind
OK_STATUSES = (STATUS_PASS, STATUS_CONSTRUCTION)

# spec constraint key -> canonical property name (regexes are strings, not numeric spec values)
_CANON = [
    (r"^density_15C_kg_m3$", "density"),
    (r"^flash_point_C$", "flash"),
    (r"^freezing_point_C$", "freezing"),
    (r"^cloud_point_C$", "cloud_point"),
    (r"^pour_point_C$", "pour_point"),
    (r"^kin_visc_(m?)(\d+)C_cSt$", "viscosity"),
    (r"^nhoc_MJ_kg$", "nhoc"),
    (r"^aromatics_wt_pct$", "aromatics_wt"),
    (r"^aromatics_vol_pct$", "aromatics_vol"),
    (r"^smoke_point_mm$", "smoke_point"),
    (r"^dcn$", "dcn"),
    (r"^ron$", "ron"),
    (r"^mon$", "mon"),
    (r"^dvpe_kPa$", "dvpe"),
    (r"^distillation_T50_C$", "distillation_T50"),
    (r"^distillation_E70_pct$", "distillation_E70"),
    (r"^ethanol_vol_pct$", "ethanol_vol"),
]


class ModelUnavailable(Exception):
    """A model cannot produce a value for this pool (unverified/missing inputs, missing registry entry...)."""


@dataclass
class PropertyModel:
    name: str
    fn: Callable                                   # fn(pool: list[str], w: np.ndarray, T_C: float | None) -> float
    margin_keys: Sequence[str] = ()                # keys of uncertainty_v4.yaml margins; combined with max()
    T_valid_C: Optional[tuple] = None              # (lo, hi) temperature range the margin was validated on
    by_construction: bool = False                  # quantity defined by the recipe itself (no model error)
    validated_classes: Optional[set] = None        # classes present in this model's validation data
    margin_missing_reason: str = ""                # why margin_keys is empty (surfaced as the reason)
    notes: str = ""
    certified_kinds: Optional[set] = None          # fuel kinds (jet/diesel/gasoline) the margin was validated on; None = any
    # {"over": [...], "under": [...]}: the LOWER limit uses the "over" margin, the UPPER limit the "under" margin;
    # when set, ``margin_keys`` is ignored for margins.
    margin_keys_by_side: Optional[dict] = None
    # margin_keys_fn(T_C of the constraint) -> keys; when set it replaces ``margin_keys`` (temperature-binned margins)
    margin_keys_fn: Optional[Callable] = None


@dataclass
class ConstraintResult:
    key: str
    prop: str
    lo: Optional[float]
    hi: Optional[float]
    model: str
    prediction: float
    applied_margin: Optional[float]
    margin_keys: tuple
    status: str
    reason: str
    slack: float                                   # min normalised margin-adjusted slack (nan if not computable)
    spec_verification: str = ""
    test_method: str = ""


@dataclass
class DesignResult:
    weights: dict
    predictions: dict
    margins: dict
    passed: bool
    extrapolating: bool
    binding_constraints: list
    constraints: list = field(default_factory=list)
    spec_name: str = ""
    pool: list = field(default_factory=list)
    absent_class_fraction: dict = field(default_factory=dict)
    optimiser: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------------------------- configuration
def load_spec(spec_name: str, configs_dir: Optional[Path] = None) -> dict:
    """``spec_name`` = ``"<standard>:<family>"`` e.g. ``"astm:jet_a"``, ``"bis:gasoline_e10"``."""
    std, fam = spec_name.split(":")
    cdir = CONFIGS if configs_dir is None else Path(configs_dir)
    with open(cdir / f"specs_{std}.yaml", "r", encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    spec = dict(doc["families"][fam])
    spec["_family"] = fam
    return spec


def fuel_kind(spec: dict) -> str:
    """'jet' | 'diesel' | 'gasoline' from the spec family name (jet_a, diesel_no2, gasoline_e10, ...)."""
    return str(spec.get("_family", "")).split("_")[0]


def load_margins(path: Optional[Path] = None) -> dict:
    with open(UNCERTAINTY_YAML if path is None else path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)["margins"]


def read_applied_margin(entry) -> tuple:
    """-> (margin, unit, direction, reason).  margin is None when the entry carries no usable ``applied_margin``.

    Only ``applied_margin`` is read.  ``value`` (the raw conformal quantile, possibly the string 'inf') is ignored.
    """
    if entry is None:
        return None, None, None, "no such entry in uncertainty_v4.yaml"
    if "applied_margin" not in entry or entry["applied_margin"] is None:
        return None, None, None, "no valid margin (applied_margin absent: too few independent groups, see yaml entry)"
    m = entry["applied_margin"]
    if isinstance(m, str) or not math.isfinite(float(m)):
        return None, None, None, "no valid margin (applied_margin is not a finite number)"
    return float(m), entry.get("unit"), entry.get("direction"), ""


def canonical_property(key: str) -> tuple:
    """-> (canonical name or None, temperature in C or None)."""
    for pat, name in _CANON:
        mo = re.match(pat, key)
        if mo:
            if name == "viscosity":
                T = float(mo.group(2))
                return name, (-T if mo.group(1) else T)
            return name, None
    return None, None


# ----------------------------------------------------------------------------------------------- evaluation helpers
def _adjust(pred: float, m: float, ln_space: bool, sign: int) -> float:
    """Worst-case value: sign=+1 -> pred + margin (upper), sign=-1 -> pred - margin (lower)."""
    if ln_space:
        return pred * math.exp(sign * m)
    return pred + sign * m


def _scale(lo, hi) -> float:
    if lo is not None and hi is not None and hi != lo:
        return abs(hi - lo)
    b = lo if lo is not None else hi
    return abs(b) if b else 1.0


def viscosity_t_bin(T_C: float) -> str:
    """Temperature bin of the lab raw-Refutas viscosity margins (same bins as scripts/v4/eval_lab_known_composition.t_bin).

    -10 C (263.15 K) and 10 C (283.15 K) belong to the WARMER bin.  The cut-offs are -10.0 C and 10.0 C
    (bins T_lt_263K | T_263_283K | T_ge_283K).
    """
    if T_C < -10.0:
        return "T_lt_263K"
    if T_C < 10.0:
        return "T_263_283K"
    return "T_ge_283K"


def _keys_for(pm: PropertyModel, T_C, side: Optional[str]):
    """-> (keys list, reason).  ``side`` is 'over' / 'under' / None."""
    if pm.margin_keys_by_side is not None:
        keys = list(pm.margin_keys_by_side.get(side, [])) if side else []
        return keys, pm.margin_missing_reason or "no uncertainty_v4.yaml margin exists for this side of the limit"
    if pm.margin_keys_fn is not None:
        if T_C is None:
            return [], "margin is temperature-binned but the constraint has no temperature"
        return list(pm.margin_keys_fn(T_C)), pm.margin_missing_reason
    return list(pm.margin_keys), pm.margin_missing_reason


def _resolved_keys(pm: PropertyModel, T_C: Optional[float]) -> tuple:
    """The yaml keys this model is judged by at ``T_C`` (for reporting)."""
    if pm.margin_keys_by_side is not None:
        return tuple(pm.margin_keys_by_side.get("over", [])) + tuple(pm.margin_keys_by_side.get("under", []))
    return tuple(_keys_for(pm, T_C, None)[0])


def _combined_margin(pm: PropertyModel, margins: dict, T_C: Optional[float] = None, side: Optional[str] = None) -> tuple:
    """-> (margin, ln_space, directions set, reason).  Several keys are combined by max (flash rule)."""
    keys, why_empty = _keys_for(pm, T_C, side)
    if not keys:
        return None, False, set(), why_empty or "no uncertainty_v4.yaml margin exists for this model"
    vals, units, dirs = [], set(), set()
    for k in keys:
        m, unit, direction, why = read_applied_margin(margins.get(k))
        if m is None:
            return None, False, set(), f"{k}: {why}"
        vals.append(m)
        units.add(unit)
        dirs.add(direction)
    return max(vals), ("ln-ratio" in units), dirs, ""


def _call(pm: PropertyModel, pool, w, T_C):
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        v = pm.fn(list(pool), np.asarray(w, float), T_C)
    ood = [str(x.message) for x in rec if x.category.__name__.endswith("OutOfDomainWarning")]
    return float(v), ood


def absent_fraction(pool, w, classes: dict, validated) -> float:
    return float(sum(wi for n, wi in zip(pool, w) if classes.get(n, "unknown") not in validated))


def forced_extrapolation(pool, classes: dict, validated) -> bool:
    """True when NO compound of the pool belongs to a validated class (the absent mass cannot be avoided)."""
    return not any(classes.get(n, "unknown") in validated for n in pool)


def evaluate_constraint(key, c: dict, pool, w, models: dict, margins: dict, pm_override: Optional[PropertyModel] = None,
                        kind: Optional[str] = None) -> ConstraintResult:
    """Evaluate one spec constraint for blend ``w`` with the primary (or an override) property model."""
    lo, hi = c.get("min"), c.get("max")
    prop, T_C = canonical_property(key)
    ver = str(c.get("verification", ""))
    meth = str(c.get("test_method", ""))
    base = dict(key=key, prop=prop or "?", lo=lo, hi=hi, spec_verification=ver, test_method=meth)
    pm = pm_override if pm_override is not None else models["props"].get(prop)
    if pm is None:
        return ConstraintResult(model="none", prediction=float("nan"), applied_margin=None, margin_keys=(),
                                status=STATUS_UNCERTIFIED, slack=float("nan"),
                                reason="no v4 model for this property", **base)
    try:
        pred, ood = _call(pm, pool, w, T_C)
    except ModelUnavailable as e:
        return ConstraintResult(model=pm.name, prediction=float("nan"), applied_margin=None,
                                margin_keys=_resolved_keys(pm, T_C), status=STATUS_UNCERTIFIED, slack=float("nan"),
                                reason=f"model unavailable: {e}", **base)
    if not math.isfinite(pred):
        return ConstraintResult(model=pm.name, prediction=pred, applied_margin=None, margin_keys=_resolved_keys(pm, T_C),
                                status=STATUS_UNCERTIFIED, slack=float("nan"),
                                reason="non-finite prediction", **base)

    # per-side margins: lower limit <- "over" margin, upper limit <- "under" margin (same margin on both sides unless
    # the model declares margin_keys_by_side)
    sided = pm.margin_keys_by_side is not None
    msides = {}
    for bound, sd in ((lo, "lo"), (hi, "hi")):
        if bound is None:
            continue
        if pm.by_construction:
            msides[sd] = (0.0, False, set(), "")
        else:
            msides[sd] = _combined_margin(pm, margins, T_C, (("over" if sd == "lo" else "under") if sided else None))
    finite = [v[0] for v in msides.values() if v[0] is not None]
    m = None if (any(v[0] is None for v in msides.values()) or not finite) else max(finite)
    sc = _scale(lo, hi)
    slacks = []
    if lo is not None:
        ml = msides["lo"]
        slacks.append((_adjust(pred, 0.0 if ml[0] is None else ml[0], ml[1], -1) - lo) / sc)
    if hi is not None:
        mh = msides["hi"]
        slacks.append((hi - _adjust(pred, 0.0 if mh[0] is None else mh[0], mh[1], +1)) / sc)
    slack = min(slacks) if slacks else float("nan")
    nominal = []
    if lo is not None:
        nominal.append(pred >= lo)
    if hi is not None:
        nominal.append(pred <= hi)
    nominal_ok = all(nominal)

    reasons = []
    if ood:
        reasons.append("model emitted OutOfDomainWarning: " + "; ".join(sorted(set(ood))))
    # extrapolation relative to this model's validation data
    validated = pm.validated_classes if pm.validated_classes is not None else models.get("validated_classes")
    if validated is not None and not pm.by_construction:
        classes = models.get("classes", {})
        af = absent_fraction(pool, w, classes, validated)
        if af > EXTRAP_WT_THRESHOLD or forced_extrapolation(pool, classes, validated):
            reasons.append(f"extrapolating: {af:.2f} wt fraction in classes absent from this model's validation data")
    if pm.T_valid_C is not None and T_C is not None and not (pm.T_valid_C[0] <= T_C <= pm.T_valid_C[1]):
        reasons.append(f"T={T_C:g} C outside the margin's validated range {pm.T_valid_C}")
    if not pm.by_construction:
        if m is None:
            for sd in ("lo", "hi"):
                if sd in msides:
                    m_sd, _ln, _dirs, why_sd = msides[sd]
                    if m_sd is None:
                        reasons.append(why_sd)
            reasons = list(dict.fromkeys(reasons))
        else:
            bad = []
            for sd, n in (("lo", "over"), ("hi", "under")):   # a minimum limit is endangered by over-prediction, a maximum by under-
                if sd in msides and not all(d in (n, "two_sided") for d in msides[sd][2]):
                    bad.append((n, sorted(map(str, msides[sd][2]))))
            if bad:
                reasons.append("margin direction does not protect the side of the limit: "
                               + "; ".join(f"{n} side has {d}" for n, d in bad))
    if (ver and ver.strip().lower() != "verified") or c.get("unverified") is True:
        reasons.append(f"spec limit verification: {ver or 'unverified: true'} (an unverified limit is never certified)")

    blocking = list(reasons)
    if pm.by_construction:
        status = STATUS_CONSTRUCTION if nominal_ok else STATUS_FAIL
        if not nominal_ok:
            reasons.append("recipe violates the limit")
        elif blocking:                                   # unverified spec limit: never certified, even by construction
            status = STATUS_UNCERTIFIED
    elif not nominal_ok or (slack == slack and slack < 0 and m is not None):
        status = STATUS_FAIL
        reasons.append("limit violated" + ("" if m is None else " after applying the margin"))
    elif blocking:
        status = STATUS_UNCERTIFIED
    else:
        status = STATUS_PASS
        if pm.certified_kinds is not None and kind not in pm.certified_kinds:
            status = STATUS_SCOPE
            reasons.append(f"margin validated on {sorted(pm.certified_kinds)} fuels only; this is a {kind or 'unknown'} spec, "
                           "so the pass is not a certification")
    return ConstraintResult(model=pm.name, prediction=pred, applied_margin=m, margin_keys=_resolved_keys(pm, T_C),
                            status=status, reason="; ".join(reasons), slack=float(slack), **base)


def evaluate_blend(spec: dict, pool, w, models: dict, margins: dict) -> list:
    kind = fuel_kind(spec)
    return [evaluate_constraint(k, c, pool, w, models, margins, kind=kind) for k, c in spec["constraints"].items()]


# ----------------------------------------------------------------------------------------------- optimiser
class _Sides:
    """Normalised margin-adjusted slack of every optimisable constraint side, cached per weight vector."""

    def __init__(self, spec, pool, models, margins):
        self.pool, self.models, self.margins, self.spec = list(pool), models, margins, spec
        self.items = []                      # (key, c, pm, T_C, (m, ln) lower side, (m, ln) upper side, sc)
        self.excluded = {}
        for key, c in spec["constraints"].items():
            prop, T_C = canonical_property(key)
            pm = models["props"].get(prop)
            if pm is None:
                self.excluded[key] = "no model"
                continue
            try:
                pm.fn(self.pool, np.full(len(self.pool), 1.0 / len(self.pool)), T_C)
            except ModelUnavailable as e:
                self.excluded[key] = f"unavailable: {e}"
                continue
            except Exception as e:           # a wrapper bug must not silently drop a constraint
                self.excluded[key] = f"error: {type(e).__name__}: {e}"
                continue
            sided = pm.margin_keys_by_side is not None
            ms = {}
            for sd, nm in (("lo", "over"), ("hi", "under")):
                r = (0.0, False, set(), "") if pm.by_construction else _combined_margin(pm, margins, T_C, nm if sided else None)
                ms[sd] = (0.0 if r[0] is None else r[0], r[1])
            self.items.append((key, c, pm, T_C, ms["lo"], ms["hi"], _scale(c.get("min"), c.get("max"))))
        self.extrap_sets = []                # distinct validated sets of the models in use (not forced)
        seen = []
        for key, c, pm, *_ in self.items:
            v = pm.validated_classes if pm.validated_classes is not None else models.get("validated_classes")
            if v is None or pm.by_construction:
                continue
            fv = frozenset(v)
            if fv in seen:
                continue
            seen.append(fv)
            if not forced_extrapolation(self.pool, models.get("classes", {}), fv):
                self.extrap_sets.append(fv)
        self._k, self._v = None, None

    def vector(self, w) -> np.ndarray:
        w = np.asarray(w, float)
        kb = w.tobytes()
        if kb == self._k:
            return self._v
        out = []
        for key, c, pm, T_C, (m_lo, ln_lo), (m_hi, ln_hi), sc in self.items:
            try:
                p = float(pm.fn(self.pool, w, T_C))
            except Exception:
                p = float("nan")
            lo, hi = c.get("min"), c.get("max")
            for bound, sgn in ((lo, -1), (hi, +1)):
                if bound is None:
                    continue
                if not math.isfinite(p):
                    out.append(BAD_SLACK)
                elif sgn < 0:
                    out.append((_adjust(p, m_lo, ln_lo, -1) - bound) / sc)
                else:
                    out.append((bound - _adjust(p, m_hi, ln_hi, +1)) / sc)
        classes = self.models.get("classes", {})
        for v in self.extrap_sets:
            out.append((EXTRAP_WT_THRESHOLD - absent_fraction(self.pool, w, classes, v)) / EXTRAP_WT_THRESHOLD)
        self._k, self._v = kb, np.array(out, float)
        return self._v


def _optimise(sides: _Sides, n: int, n_starts: int, seed: int) -> tuple:
    rng = np.random.default_rng(seed)
    starts = [np.full(n, 1.0 / n)] + [rng.dirichlet(np.ones(n)) for _ in range(max(n_starts - 1, 0))]
    n_sides = len(sides.vector(starts[0]))
    best_w, best_score, n_ok = starts[0], -np.inf, 0
    if n_sides == 0:
        return best_w, float("nan"), 0
    for w0 in starts:
        s0 = float(np.min(sides.vector(w0)))
        z0 = np.append(w0, s0)
        cons = [{"type": "eq", "fun": lambda z: np.sum(z[:-1]) - 1.0},
                {"type": "ineq", "fun": lambda z: sides.vector(np.clip(z[:-1], 0, None)) - z[-1]}]
        bnds = [(0, 1)] * n + [(None, None)]
        try:
            res = minimize(lambda z: -z[-1], z0, method="SLSQP", bounds=bnds, constraints=cons,
                           options={"maxiter": MAX_ITER, "ftol": 1e-9})
            w = np.clip(res.x[:-1], 0, None)
            w = w / w.sum()
            n_ok += int(res.success)
        except Exception:
            w = w0
        for cand in (w, w0):
            sc = float(np.min(sides.vector(cand)))
            if sc > best_score:
                best_score, best_w = sc, cand
    return best_w, best_score, n_ok


def design(spec_name: str, pool: list, models: dict, n_starts: int = 30, seed: int = 42,
           margins_path: Optional[Path] = None, configs_dir: Optional[Path] = None) -> DesignResult:
    """Search the simplex of ``pool`` for the blend with the largest worst-case normalised slack, then report the
    honest status of every constraint (see module docstring)."""
    spec = load_spec(spec_name, configs_dir)
    margins = load_margins(margins_path)
    pool = list(pool)
    sides = _Sides(spec, pool, models, margins)
    w, score, n_ok = _optimise(sides, len(pool), n_starts, seed)
    w = np.where(w < 1e-6, 0.0, w)
    w = w / w.sum()

    cons = evaluate_blend(spec, pool, w, models, margins)
    for cr in cons:
        if cr.key in sides.excluded and not cr.reason:
            cr.reason = sides.excluded[cr.key]

    classes = models.get("classes", {})
    absent, extrap = {}, False
    used = {}
    for cr in cons:
        pm = models["props"].get(cr.prop)
        if pm is None or pm.by_construction:
            continue
        v = pm.validated_classes if pm.validated_classes is not None else models.get("validated_classes")
        if v is None:
            continue
        af = absent_fraction(pool, w, classes, v)
        used[pm.name] = af
        if af > EXTRAP_WT_THRESHOLD or forced_extrapolation(pool, classes, v):
            extrap = True
    absent = used

    all_ok = all(cr.status in OK_STATUSES for cr in cons)
    passed = bool(all_ok and not extrap)

    binding = []
    for cr in cons:
        if cr.status not in OK_STATUSES:
            binding.append(dict(constraint=cr.key, status=cr.status, reason=cr.reason, slack=cr.slack))
    if extrap:
        binding.append(dict(constraint="extrapolation", status=STATUS_FAIL,
                            reason=f"> {EXTRAP_WT_THRESHOLD} wt fraction in classes absent from validation data: {absent}",
                            slack=float("nan")))
    if not binding and cons:
        finite = [cr for cr in cons if cr.slack == cr.slack]
        if finite:
            tight = min(finite, key=lambda cr: cr.slack)
            binding.append(dict(constraint=tight.key, status="tightest_passing", reason=tight.reason, slack=tight.slack))

    return DesignResult(
        weights={n: float(x) for n, x in zip(pool, w)},
        predictions={cr.key: cr.prediction for cr in cons},
        margins={cr.key: cr.applied_margin for cr in cons},
        passed=passed, extrapolating=bool(extrap), binding_constraints=binding, constraints=cons,
        spec_name=spec_name, pool=pool, absent_class_fraction=absent,
        optimiser=dict(best_min_slack=float(score), n_starts=int(n_starts), n_converged=int(n_ok), seed=int(seed),
                       excluded_from_optimiser=dict(sides.excluded)))


# ----------------------------------------------------------------------------------------------- robustness
def robustness_table(spec_name: str, pool: list, weights, models: dict, margins_path: Optional[Path] = None,
                     configs_dir: Optional[Path] = None) -> list:
    """Evaluate every constraint with EVERY available model form (primary + ``models['alternates'][prop]``).

    Each model function is actually called (once per constraint); nothing is added to a prediction as an offset.
    Each model form is judged with its OWN yaml margin (or reported uncertified when it has none)."""
    spec = load_spec(spec_name, configs_dir)
    margins = load_margins(margins_path)
    w = np.asarray(weights, float)
    rows = []
    for key, c in spec["constraints"].items():
        prop, _ = canonical_property(key)
        forms = []
        if prop in models["props"]:
            forms.append((True, models["props"][prop]))
        forms += [(False, pm) for pm in models.get("alternates", {}).get(prop, [])]
        for primary, pm in forms:
            cr = evaluate_constraint(key, c, pool, w, models, margins, pm_override=pm, kind=fuel_kind(spec))
            rows.append(dict(spec=spec_name, constraint=key, model=pm.name, is_primary=primary,
                             prediction=cr.prediction, applied_margin=cr.applied_margin, lo=cr.lo, hi=cr.hi,
                             status=cr.status, reason=cr.reason))
    return rows
