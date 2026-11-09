"""Provenance ledger for the real-fuel tables (v4 Task 1).

The v2/v3 real-fuel tables (``data/real_fuels/inventory*.csv``,
``compositions*.csv`` and ``data/real_fuels/v3/*``) mix three kinds of numbers:
values transcribed from a source file, values typed into an ingest script, and
values *computed* by the ingest script (from formulas, from the fuel's own
target properties, or as one constant pasted onto many fuels). Their labels
(``composition_type="measured_dha_gcxgc"``, ``source=...``) do not tell these
apart, so this module never trusts them. It classifies every value from
evidence:

* static evidence -- the ingest scripts are parsed with :mod:`ast` and the
  code that builds each table row is *replayed* in a restricted evaluator
  (literal fuel tables, arithmetic, ``min/max/round``; no I/O, no attribute
  access except ``dict.get`` and ``np.nan``). Each value gets the expression
  that produced it and the set of per-fuel input variables it depends on.
* data evidence -- exact functional relations tested on the table values
  (``freeze == cloud - 3``, FACE ``toluene == clip((dens-680)/200)`` after the
  script's normalisation, identical fraction vectors across fuels, constants
  shared by many fuels) and per-fuel, per-property value matches against the
  raw files under ``extras/fuelblend_data/data/raw`` and ``data/raw``.

Classes (``CLASSES``):

``raw_file``
    value equals the value for *this fuel and this property* in a raw file
    (never a global "the number appears somewhere" search).
``cited_typed``
    a per-fuel literal typed into the script (or read from a processed file
    that carries a per-row reference), with no formula, constant-sharing or
    raw-file conflict detected. Not verified against a raw file.
``formula_derived``
    computed by a formula from other values (including model reconstructions;
    their evidence starts with ``reconstructed:``).
``target_derived_composition``
    a composition fraction computed from the same fuel's own property values
    (RON, MON, density, DCN, ...) -- circular for any property model.
``placeholder_constant``
    one literal shared by two or more fuels (or one fraction vector shared by
    two or more systems), not backed by the raw file for that fuel.
``unknown``
    none of the above could be shown (including typed values that *conflict*
    with the raw file the script claims as source, and structural zeros).
"""
from __future__ import annotations

import ast
import math
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

CLASSES = ("raw_file", "cited_typed", "formula_derived", "target_derived_composition",
           "placeholder_constant", "unknown")
LEDGER_COLUMNS = ["fuel_id", "field", "value", "provenance_class", "evidence"]

INGEST_SCRIPTS = ("scripts/v2/ingest_real_fuels_v2.py", "scripts/v3/ingest_phase1_expansion.py")
DEFAULT_RAW_ROOTS = ("extras/fuelblend_data/data/raw", "data/raw")

PROPERTY_FIELDS = ("ron", "mon", "dvpe_kpa", "density_kg_m3", "viscosity_kin_cSt", "flash_point_C",
                   "freezing_point_C", "cloud_point_C", "pour_point_C", "dcn", "heat_of_combustion_MJ_kg")
META_FIELDS = ("number_of_resolved_components", "unresolved_fraction",
               "ethanol_vol_pct", "ethanol_mass_pct", "oxygenate_pct")
STRUCTURAL_ZERO_FIELDS = ("unresolved_fraction", "ethanol_vol_pct", "ethanol_mass_pct", "oxygenate_pct")
NUMERIC_FIELDS = PROPERTY_FIELDS + META_FIELDS
# Loop variables of the ingest scripts' per-fuel literal tables that carry a target property.
# (Names as written in scripts/v2/ingest_real_fuels_v2.py and scripts/v3/ingest_phase1_expansion.py;
# in v3 the same variables are also written to the property columns, e.g. "ron": ron.)
PROPERTY_VAR_NAMES = {"ron": "ron", "mon": "mon", "ron_val": "ron", "mon_val": "mon",
                      "dens": "density_kg_m3", "dens_val": "density_kg_m3", "dvpe": "dvpe_kpa",
                      "dvpe_val": "dvpe_kpa", "dcn": "dcn", "visc40": "viscosity_kin_cSt",
                      "visc": "viscosity_kin_cSt", "flash": "flash_point_C", "freeze": "freezing_point_C",
                      "cloud": "cloud_point_C", "pour": "pour_point_C"}

_REPO_ROOT = Path(__file__).resolve().parents[3]
_UNK = object()  # "could not evaluate"


# =========================================================================== small helpers
def _fmt(v: Any) -> str:
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return format(float(v), ".10g")
    return str(v)


def _num(v: Any) -> float | None:
    if isinstance(v, (bool, np.bool_)):
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _close(a: Any, b: Any, rel: float = 1e-9, abs_: float = 1e-9) -> bool:
    fa, fb = _num(a), _num(b)
    if fa is None or fb is None:
        return False
    return abs(fa - fb) <= max(abs_, rel * max(abs(fa), abs(fb)))


def _same(a: Any, b: Any) -> bool:
    if _num(a) is not None and _num(b) is not None:
        return _close(a, b, rel=1e-9, abs_=5e-7)
    if isinstance(a, (str, bool, np.bool_)) or isinstance(b, (str, bool, np.bool_)):
        return str(a) == str(b)
    return False


def _relpath(p: Path, root: Path) -> str:
    try:
        return p.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return p.as_posix()


# =========================================================================== static replay
class _NanStub:
    nan = float("nan")


class _DefaultGet:
    """Stands in for lookup dicts built from files (e.g. ``smi2ik``): ``.get(k, d) -> d``."""

    def get(self, _k, default=None):
        return default


_SAFE_BUILTINS = {"min": min, "max": max, "round": round, "sum": sum, "abs": abs, "int": int,
                  "float": float, "len": len, "str": str, "bool": bool, "zip": zip, "range": range}
_ALLOWED_NODES = (ast.Expression, ast.Constant, ast.Name, ast.Load, ast.Store, ast.BinOp, ast.UnaryOp,
                  ast.BoolOp, ast.Compare, ast.IfExp, ast.Call, ast.Attribute, ast.Subscript, ast.Slice,
                  ast.Tuple, ast.List, ast.Dict, ast.ListComp, ast.comprehension, ast.keyword,
                  ast.operator, ast.unaryop, ast.boolop, ast.cmpop, ast.JoinedStr, ast.FormattedValue)


_SAFE_ATTRS = ("get", "nan", "replace", "strip", "lower", "upper")


def _safe_eval(node: ast.AST, env: dict) -> Any:
    for sub in ast.walk(node):
        if not isinstance(sub, _ALLOWED_NODES):
            return _UNK
        if isinstance(sub, ast.Attribute) and sub.attr not in _SAFE_ATTRS:
            return _UNK
        if isinstance(sub, ast.Call):
            f = sub.func
            ok = (isinstance(f, ast.Name) and f.id in _SAFE_BUILTINS) or (
                isinstance(f, ast.Attribute) and f.attr in _SAFE_ATTRS and f.attr != "nan")
            if not ok:
                return _UNK
    try:
        code = compile(ast.Expression(body=node), "<replay>", "eval")
        val = eval(code, {"__builtins__": _SAFE_BUILTINS}, env)  # noqa: S307 - whitelisted AST only
    except Exception:  # noqa: BLE001 - any failure means "not evaluable"
        return _UNK
    return val


def _names(node: ast.AST) -> set[str]:
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def _is_const(node: ast.AST) -> bool:
    if isinstance(node, ast.Constant):
        return True
    if isinstance(node, ast.UnaryOp) and isinstance(node.operand, ast.Constant):
        return True
    return False


def _ifexp_leaf(node: ast.AST, env: dict) -> ast.AST | None:
    """Follow an IfExp chain of constants with the current env; return the leaf taken."""
    while isinstance(node, ast.IfExp):
        t = _safe_eval(node.test, env)
        if t is _UNK:
            return None
        node = node.body if t else node.orelse
    return node if _is_const(node) else None


def _is_const_ifexp(node: ast.AST) -> bool:
    if not isinstance(node, ast.IfExp):
        return False
    stack = [node]
    while stack:
        n = stack.pop()
        if isinstance(n, ast.IfExp):
            stack += [n.body, n.orelse]
        elif not _is_const(n):
            return False
    return True


@dataclass
class CodeRecord:
    """One dict literal appended to an output list, as replayed for one fuel (or as a template)."""

    script: str
    loop_line: int            # line of the innermost enclosing ``for`` (0 if none)
    outer_loop_line: int      # line of the outermost enclosing ``for``
    line: int
    kind: str                 # 'inventory' | 'composition'
    replayed: bool            # every enclosing loop iterated literal data
    keys: tuple[str, ...]
    nodes: dict[str, ast.AST]
    values: dict[str, Any]
    segs: dict[str, str]
    deps: dict[str, frozenset]
    env: dict                 # snapshot for late (per-fuel) evaluation of templates
    loop_vars: frozenset
    carry_origin: dict[str, str] = field(default_factory=dict)   # key -> table path read by the loop
    model_calls: tuple[str, ...] = ()                            # non-whitelisted calls in the loop
    ordinal: int = -1

    @property
    def fuel_id(self) -> str | None:
        for k in ("fuel_id", "system_id"):
            v = self.values.get(k, _UNK)
            if isinstance(v, str):
                return v[5:] if v.startswith("fuel:") else v
        return None

    @property
    def template_id(self) -> str:
        return f"{self.script}:L{self.line}"


class _Replayer:
    """Replays the row-building code of one ingest script without running any I/O."""

    def __init__(self, script: str, code: str):
        self.script = script
        self.code = code
        self.tree = ast.parse(code)
        self.records: list[CodeRecord] = []
        self.defs: dict[str, tuple[str, int, frozenset]] = {}   # dep-tag -> (segment, line, direct deps)
        self.prop_tags: dict[str, str] = {}                      # "var@line" -> property key
        self.literal_defs: set[str] = set()                      # dep-tags assigned a pure literal
        self.paths = self._module_paths()
        self._index_property_vars()

    # ------------------------------------------------------------------ static indexes
    def _path_of(self, n: ast.AST, known: dict[str, str] | None = None) -> str | None:
        """ROOT / "a" / "b" (or KNOWN / "c") -> "a/b" ; None if not a constant path expression."""
        known = self.paths if known is None else known

        def parts(x: ast.AST) -> list[str] | None:
            if isinstance(x, ast.BinOp) and isinstance(x.op, ast.Div):
                left, right = parts(x.left), parts(x.right)
                return None if left is None or right is None else left + right
            if isinstance(x, ast.Name):
                if x.id == "ROOT":
                    return []
                return known[x.id].split("/") if x.id in known else None
            if isinstance(x, ast.Constant) and isinstance(x.value, str):
                return [x.value]
            return None

        p = parts(n)
        return "/".join(s for s in p if s) if p else None

    def _module_paths(self) -> dict[str, str]:
        """Module-level NAME = ROOT / "a" / "b"  ->  {"NAME": "a/b"}."""
        out: dict[str, str] = {}
        for st in self.tree.body:
            if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], ast.Name):
                p = self._path_of(st.value, out)
                if p:
                    out[st.targets[0].id] = p
        return out

    def _index_property_vars(self) -> None:
        """Per-fuel loop variables that hold a fuel *property* (by the script's own variable name)."""
        for node in ast.walk(self.tree):
            if not isinstance(node, ast.For):
                continue
            for n in _names(node.target):
                if n in PROPERTY_VAR_NAMES:
                    self.prop_tags[f"{n}@{node.lineno}"] = PROPERTY_VAR_NAMES[n]

    def seg(self, node: ast.AST) -> str:
        s = ast.get_source_segment(self.code, node) or ast.unparse(node)
        return " ".join(s.split())

    # ------------------------------------------------------------------ replay
    def run(self) -> "_Replayer":
        env: dict = {"np": _NanStub(), "smi2ik": _DefaultGet(), "name2ik": _DefaultGet()}
        for fn in self.tree.body:
            if isinstance(fn, ast.FunctionDef) and fn.name == "main":
                self._exec(fn.body, env, {}, {}, [], True)
        return self

    def _deps(self, node: ast.AST, deps: dict) -> frozenset:
        out: set = set()
        for n in _names(node):
            out |= deps.get(n, frozenset())
        return frozenset(out)

    def _origin(self, node: ast.AST, origins: dict) -> str | None:
        """Table a value is read from: follows the *data* operand (subscript/attribute base, call
        argument, arithmetic operand), never index expressions such as ``df[df.id == sys_id]``."""
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "read_csv" and node.args and isinstance(node.args[0], ast.Name)):
            return self.paths.get(node.args[0].id)
        if isinstance(node, ast.Name):
            return origins.get(node.id)
        if isinstance(node, (ast.Subscript, ast.Attribute)):
            return self._origin(node.value, origins)
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute):
                return self._origin(node.func.value, origins)
            if isinstance(node.func, ast.Name) and node.args:
                for a in node.args:
                    o = self._origin(a, origins)
                    if o:
                        return o
            return None
        if isinstance(node, ast.BinOp):
            return self._origin(node.left, origins) or self._origin(node.right, origins)
        return None

    def _exec(self, stmts, env, deps, origins, loops, literal) -> None:
        for st in stmts:
            if isinstance(st, ast.Assign) and len(st.targets) == 1:
                self._assign(st, st.targets[0], env, deps, origins)
            elif isinstance(st, ast.For):
                self._for(st, env, deps, origins, loops, literal)
            elif isinstance(st, ast.If):
                t = _safe_eval(st.test, env)
                if t is _UNK:
                    self._exec(st.body, env, deps, origins, loops, literal)
                    self._exec(st.orelse, env, deps, origins, loops, literal)
                    # either branch may have run: anything assigned under the If is unknown afterwards
                    for sub in ast.walk(st):
                        if isinstance(sub, ast.Assign):
                            for tg in sub.targets:
                                if isinstance(tg, ast.Name):
                                    env[tg.id] = _UNK
                else:
                    self._exec(st.body if t else st.orelse, env, deps, origins, loops, literal)
            elif isinstance(st, (ast.With, ast.Try)):
                self._exec(st.body, env, deps, origins, loops, literal)
            elif isinstance(st, ast.Expr) and isinstance(st.value, ast.Call):
                c = st.value
                if (isinstance(c.func, ast.Attribute) and c.func.attr == "append"
                        and len(c.args) == 1 and isinstance(c.args[0], ast.Dict)):
                    self._record(c.args[0], env, deps, origins, loops, literal)

    def _assign(self, st, tgt, env, deps, origins) -> None:
        val = _safe_eval(st.value, env)
        d = self._deps(st.value, deps)
        org = self._origin(st.value, origins)
        if isinstance(tgt, ast.Name) and isinstance(st.value, ast.BinOp) and isinstance(st.value.op, ast.Div):
            p = self._path_of(st.value)
            if p:
                self.paths[tgt.id] = p
        if isinstance(tgt, ast.Tuple):
            for n in _names(tgt):
                env[n] = _UNK
                deps[n] = d | {f"{n}@{st.lineno}"}
                self.defs[f"{n}@{st.lineno}"] = (self.seg(st), st.lineno, d)
                if org:
                    origins[n] = org
            return
        if isinstance(tgt, ast.Name):
            tag = f"{tgt.id}@{st.lineno}"
            self.defs[tag] = (self.seg(st), st.lineno, d)
            if all(isinstance(n, (ast.Constant, ast.Name, ast.Tuple, ast.List, ast.Load, ast.UnaryOp,
                                  ast.unaryop)) for n in ast.walk(st.value)):
                self.literal_defs.add(tag)   # a literal table / constant (names only refer to constants)
            env[tgt.id] = val
            deps[tgt.id] = d | {tag}
            if org:
                origins[tgt.id] = org
            else:
                origins.pop(tgt.id, None)
        elif isinstance(tgt, ast.Subscript):
            root = tgt
            while isinstance(root, ast.Subscript):
                root = root.value
            container, key = _safe_eval(tgt.value, env), _safe_eval(tgt.slice, env)
            if container is not _UNK and key is not _UNK and val is not _UNK:
                try:
                    container[key] = val
                except Exception:  # noqa: BLE001
                    pass
            if isinstance(root, ast.Name):
                deps[root.id] = deps.get(root.id, frozenset()) | d

    def _for(self, st, env, deps, origins, loops, literal) -> None:
        it = _safe_eval(st.iter, env)
        tnames = _names(st.target)
        base = self._deps(st.iter, deps)
        org = self._origin(st.iter, origins)
        loop_info = (st.lineno, tnames)
        for n in tnames:
            deps[n] = base | {f"{n}@{st.lineno}"}
            if org:
                origins[n] = org
            else:
                origins.pop(n, None)
        if it is _UNK or not isinstance(it, (list, tuple)):
            for n in tnames:
                env[n] = _UNK
            self._exec(st.body, env, deps, origins, loops + [loop_info], False)
            return
        for item in it:
            try:
                self._bind(st.target, item, env)
            except Exception:  # noqa: BLE001
                for n in tnames:
                    env[n] = _UNK
            self._exec(st.body, env, deps, origins, loops + [loop_info], literal)

    def _bind(self, tgt, item, env) -> None:
        if isinstance(tgt, ast.Name):
            env[tgt.id] = item
        elif isinstance(tgt, ast.Tuple):
            if len(tgt.elts) != len(item):
                raise ValueError
            for t, v in zip(tgt.elts, item):
                self._bind(t, v, env)
        else:
            raise ValueError

    def _model_calls(self, loop_line: int) -> tuple[str, ...]:
        for node in ast.walk(self.tree):
            if isinstance(node, ast.For) and node.lineno == loop_line:
                calls = set()
                for sub in ast.walk(node):
                    if isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name) \
                            and sub.func.id not in _SAFE_BUILTINS:
                        calls.add(sub.func.id)
                return tuple(sorted(calls))
        return ()

    def _record(self, d: ast.Dict, env, deps, origins, loops, literal) -> None:
        keys, nodes, values, segs, kdeps, carry = [], {}, {}, {}, {}, {}
        for k, v in zip(d.keys, d.values):
            if not (isinstance(k, ast.Constant) and isinstance(k.value, str)):
                continue
            key = k.value
            keys.append(key)
            nodes[key] = v
            values[key] = _safe_eval(v, env)
            segs[key] = self.seg(v)
            kdeps[key] = self._deps(v, deps)
            org = self._origin(v, origins)
            if org and isinstance(v, (ast.Subscript, ast.Call)):
                carry[key] = org
        kind = "composition" if "inchikey" in keys else "inventory"
        loop_vars = frozenset(n for _, ns in loops for n in ns)
        loop_line = loops[-1][0] if loops else 0
        outer = loops[0][0] if loops else 0
        self.records.append(CodeRecord(
            script=self.script, loop_line=loop_line, outer_loop_line=outer, line=d.lineno, kind=kind,
            replayed=bool(loops) and literal, keys=tuple(keys), nodes=nodes, values=values, segs=segs,
            deps=kdeps, env=dict(env), loop_vars=loop_vars, carry_origin=carry,
            model_calls=self._model_calls(outer) if outer else ()))

    # ------------------------------------------------------------------ dependency explanation
    def property_deps(self, dset: frozenset) -> list[str]:
        return sorted({self.prop_tags[t] for t in dset if t in self.prop_tags})

    def explain(self, dset: frozenset, limit: int = 6) -> list[str]:
        """Assignments in the dependency chain whose own expression reads a per-fuel property var."""
        out = []
        for tag in sorted(dset, key=lambda t: (self.defs.get(t, ("", 10 ** 9))[1], t)):
            if tag not in self.defs:
                continue
            seg, line, direct = self.defs[tag]
            if any(t in self.prop_tags for t in direct):
                out.append(f"{self.script}:L{line} `{seg}`")
        return out[:limit]


@lru_cache(maxsize=8)
def _replay_cached(script: str, code: str) -> _Replayer:
    return _Replayer(script, code).run()


def load_source_code(root: Path | str = _REPO_ROOT) -> dict[str, str]:
    root = Path(root)
    return {s: (root / s).read_text(encoding="utf-8") for s in INGEST_SCRIPTS if (root / s).exists()}


def scan_source(source_code: dict[str, str]) -> list[CodeRecord]:
    """Static pass: every dict appended to an output list in every ingest script."""
    recs: list[CodeRecord] = []
    for script in sorted(source_code):
        recs.extend(_replay_cached(script, source_code[script]).records)
    return recs


# =========================================================================== raw files
class RawIndex:
    """Per-fuel, per-property lookups into the raw files (never a global value search)."""

    CRC647_FIELDS = {"density_kg_m3": "density_kg_m3", "flash_point_C": "flash_point_C",
                     "freezing_point_C": "freezing_point_C", "pour_point_C": "pour_point_C",
                     "viscosity_kin_cSt": "viscosity_kin_cSt", "distillation_T10_C": "distillation_T10_C",
                     "heat_of_combustion_net_MJ_kg": "net_heat_content_btu_lb",
                     "heat_of_combustion_MJ_kg": "net_heat_content_btu_lb"}
    COOPTIMA_NAMES = {"TIGAS_Green_Gasoline": "TIGAS Green Gasoline",
                      "CoOptima_4comp": "4 Component Surrogate", "NREL_BOB_01": "2016 NREL BOB-01"}
    BTU_LB_TO_MJ_KG = 0.002326

    def __init__(self, raw_roots, repo_root: Path):
        self.repo_root = repo_root
        self.roots = [Path(r) if Path(r).is_absolute() else repo_root / r for r in raw_roots]
        self._cache: dict[str, tuple[str, pd.DataFrame] | None] = {}

    def _find(self, rel: str) -> tuple[str, pd.DataFrame] | None:
        if rel in self._cache:
            return self._cache[rel]
        hit = None
        for r in self.roots:
            p = r / rel
            if p.exists():
                hit = (_relpath(p, self.repo_root),
                       pd.read_csv(p, dtype=str, encoding_errors="replace", keep_default_na=False))
                break
        self._cache[rel] = hit
        return hit

    # ---------------------------------------------------------------- scalar properties
    def scalar(self, fuel_id: str, prop: str) -> tuple[str, list[tuple[float, str]]] | None:
        """(file, [(value_in_ledger_units, provenance note), ...]) or None if no raw source maps here."""
        m = re.fullmatch(r"CRC647_(\d+)", fuel_id)
        if m and prop in self.CRC647_FIELDS:
            hit = self._find("pour_point/crc647_fuel_properties.csv")
            if hit is None:
                return None
            path, df = hit
            code = df["fuel_code"].str.strip().str.lstrip("0")
            sub = df[(code == m.group(1).lstrip("0")) & (df["property"] == self.CRC647_FIELDS[prop])]
            vals = []
            for _, r in sub.iterrows():
                v = _num(r["value_std"])
                if v is None:
                    continue
                if self.CRC647_FIELDS[prop] == "net_heat_content_btu_lb":
                    v = round(v * self.BTU_LB_TO_MJ_KG, 4)
                note = f"{r['property']}={r['value_std']}{' @T=' + r['T_C'] + 'C' if r['T_C'] else ''} " \
                       f"({r['method']}, {r['source_pdf']} p{r['page']})"
                vals.append((v, note))
            return path, vals
        if prop in ("ron", "mon"):
            name = None
            m = re.fullmatch(r"PRF_E(\d+)", fuel_id)
            if m:
                name, rel = f"PRF+{m.group(1)}", "octane/face_gasoline_set_I.csv"
            m = re.fullmatch(r"CRC_FACE_([A-J])_E(\d+)", fuel_id)
            if m:
                name, rel = f"FACE {m.group(1)}+{m.group(2)}", "octane/face_gasoline_set_II.csv"
            if name:
                hit = self._find(rel)
                if hit is None:
                    return None
                path, df = hit
                sub = df[df["Name"].str.strip() == name]
                col = prop.upper()
                return path, [(v, f"{name}: {col}={r[col]}") for _, r in sub.iterrows()
                              if (v := _num(r[col])) is not None]
        if fuel_id in self.COOPTIMA_NAMES and prop in ("ron", "mon", "density_kg_m3", "dvpe_kpa"):
            hit = self._find("cooptima/fuel_engine_co_optimization.csv")
            if hit is None:
                return None
            path, df = hit
            sub = df[df["Blend_Name"].str.strip() == self.COOPTIMA_NAMES[fuel_id]]
            col, scale = {"ron": ("Blend_RON", 1.0), "mon": ("Blend_MON", 1.0),
                          "density_kg_m3": ("Blend_Density", 1000.0), "dvpe_kpa": ("Blend_DVPE", 1.0)}[prop]
            return path, [(round(v * scale, 6), f"{self.COOPTIMA_NAMES[fuel_id]}: {col}={r[col]}")
                          for _, r in sub.iterrows() if (v := _num(r[col])) is not None]
        return None

    def check_scalar(self, fuel_id: str, prop: str, value: Any) -> tuple[str, str] | None:
        """('match'|'conflict'|'absent', evidence) or None when no raw source maps to (fuel, prop)."""
        got = self.scalar(fuel_id, prop)
        if got is None:
            return None
        path, vals = got
        if not vals:
            return "absent", f"raw: {path} has no {prop} for {fuel_id}"
        for v, note in vals:
            if _close(v, value, rel=1e-6, abs_=1e-6):
                return "match", f"raw: {path} -> {note}"
        listed = "; ".join(n for _, n in vals[:4])
        return "conflict", f"raw CONFLICT: {path} gives {listed}"

    # ---------------------------------------------------------------- compositions
    def composition(self, fuel_id: str) -> tuple[str, str, pd.DataFrame] | None:
        """(file, basis, df[smiles, frac]) for fuels whose raw molecule list maps to fuel_id."""
        if fuel_id == "NIST_RP1":
            hit = self._find("rp1/rp1_blends_composition.csv")
            if hit is None:
                return None
            path, df = hit
            sub = df[df["blend_id"] == "RP1_sample_1_original"]
            basis = sub["frac_type"].iloc[0] if len(sub) else "mass"
            return path, basis, pd.DataFrame({"smiles": sub["component_smiles"].str.strip(),
                                              "frac": pd.to_numeric(sub["frac_value"], errors="coerce")})
        hit = self._find("surrogates/fuel_surrogates.csv")
        if hit is not None:
            path, df = hit
            sub = df[df["fuel_code"] == fuel_id]
            if len(sub):
                return path, sub["frac_type"].iloc[0], pd.DataFrame(
                    {"smiles": sub["component_smiles"].str.strip(),
                     "frac": pd.to_numeric(sub["frac_value"], errors="coerce")})
        return None

    def gc_classes(self, fuel_id: str) -> tuple[str, pd.DataFrame] | None:
        hit = self._find("gc_composition/fuel_gc_compositions.csv")
        if hit is None:
            return None
        path, df = hit
        sub = df[df["fuel_code"] == fuel_id.replace("_", "-")]
        return (path, sub) if len(sub) else None


# =========================================================================== data relations
def _face_gasoline_fracs(ron, mon, dens, unres):
    """The v2/v3 FACE gasoline allocation, replayed from the fuel's own RON/MON/density."""
    s_val = ron - mon
    aro = min(0.45, max(0.05, (dens - 680.0) / 200.0))
    ole = min(0.25, max(0.01, s_val / 20.0))
    iso = max(0.10, (ron / 100.0) * (1.0 - aro - ole - unres))
    hep = max(0.05, 1.0 - unres - aro - ole - iso)
    raw = [iso, hep, aro, ole]
    k = (1.0 - unres) / sum(raw)
    return [round(f * k, 5) for f in raw]


def _face_diesel_fracs(dcn, dens, unres):
    aro = min(0.40, max(0.10, (dens - 800.0) / 150.0))
    npar = min(0.50, max(0.15, (dcn - 25.0) / 70.0))
    hmn = min(0.35, max(0.10, (1.0 - unres - aro - npar) * 0.6))
    dec = max(0.05, 1.0 - unres - aro - npar - hmn)
    raw = [npar, hmn, aro, dec]
    k = (1.0 - unres) / sum(raw)
    return [round(f * k, 5) for f in raw]


_GASOLINE_ORDER = ("NHTMVDHEPJAVLT-UHFFFAOYSA-N", "IMNFDUFMRHMDMM-UHFFFAOYSA-N",
                   "UAEPNZWRGJTJPN-UHFFFAOYSA-N", "YWAKXOMACGNMNL-UHFFFAOYSA-N")
_DIESEL_ORDER = ("DCAYPVUWAIABOU-UHFFFAOYSA-N", "GWHJZXXIDMPWGX-UHFFFAOYSA-N",
                 "HCAJNOARHGCSNO-UHFFFAOYSA-N", "NNBZCPXTIHJBJL-UHFFFAOYSA-N")


def _row_num(row: pd.Series, col: str) -> float | None:
    return _num(row[col]) if col in row.index else None


def _inventory_relations(inv: pd.DataFrame) -> dict[tuple[str, str], str]:
    """Exact functional relations that hold in the table values: {(fuel_id, field): description}."""
    out: dict[tuple[str, str], str] = {}
    by_id = {str(r["fuel_id"]): r for _, r in inv.iterrows()}
    for fid, r in by_id.items():
        g = lambda c: _row_num(r, c)  # noqa: E731
        fr, cl, po, fl, de = g("freezing_point_C"), g("cloud_point_C"), g("pour_point_C"), \
            g("flash_point_C"), g("density_kg_m3")
        if fr is not None and cl is not None and _close(fr, cl - 3.0):
            out[(fid, "freezing_point_C")] = f"data: freezing_point_C ({_fmt(fr)}) == cloud_point_C - 3 exactly"
        if po is not None and fr is not None and _close(po, fr - 4.0):
            out[(fid, "pour_point_C")] = f"data: pour_point_C ({_fmt(po)}) == freezing_point_C - 4 exactly"
        if fl is not None and de is not None and _close(fl, round(52.0 + (de - 800.0) * 0.25, 1)):
            out[(fid, "flash_point_C")] = (f"data: flash_point_C ({_fmt(fl)}) == round(52 + "
                                           f"(density {_fmt(de)} - 800) * 0.25, 1) exactly")
        m = re.fullmatch(r"(.+)_E(\d+)", fid)
        if m and m.group(1) in by_id:
            b = by_id[m.group(1)]
            v_eth = float(m.group(2))
            dv, bdv = g("dvpe_kpa"), _row_num(b, "dvpe_kpa")
            if dv is not None and bdv is not None:
                step = 7.0 if v_eth <= 20 else 3.0
                if _close(dv, round(bdv + step, 1)):
                    out[(fid, "dvpe_kpa")] = (f"data: dvpe_kpa ({_fmt(dv)}) == base {m.group(1)} dvpe "
                                              f"({_fmt(bdv)}) + {step} (rule: +7.0 if E<=20 else +3.0)")
            bde = _row_num(b, "density_kg_m3")
            if de is not None and bde is not None and _close(
                    de, round(((100.0 - v_eth) * bde + v_eth * 789.0) / 100.0, 1)):
                out[(fid, "density_kg_m3")] = (f"data: density_kg_m3 ({_fmt(de)}) == volume-weighted mix of base "
                                               f"{m.group(1)} ({_fmt(bde)}) and ethanol 789.0")
    return out


# =========================================================================== classifier
@dataclass
class _Ctx:
    inv: pd.DataFrame
    comp: pd.DataFrame
    recs: list[CodeRecord]
    reps: dict[str, _Replayer]
    raw: RawIndex
    relations: dict[tuple[str, str], str]
    inv_match: dict[str, CodeRecord | None] = field(default_factory=dict)
    groups: dict | None = None


def _record_value(rec: CodeRecord, key: str, fuel_id: str) -> Any:
    """Value a record assigns to ``key`` for ``fuel_id`` (templates are evaluated with fid bound)."""
    v = rec.values.get(key, _UNK)
    if v is not _UNK or rec.replayed:
        return v
    node = rec.nodes.get(key)
    if node is None:
        return _UNK
    env = dict(rec.env)
    env["fid"] = fuel_id
    if _is_const(node) or _is_const_ifexp(node):
        return _safe_eval(node, env)
    return _UNK


def _match_inventory(fid: str, row: pd.Series, columns: set[str], recs: list[CodeRecord]) -> CodeRecord | None:
    best, best_score = None, 0
    for rec in recs:
        if rec.kind != "inventory":
            continue
        if rec.replayed and rec.fuel_id != fid:
            continue
        score = 0
        for k in rec.keys:
            if k not in row.index:
                continue
            dv = row[k]
            if dv is None or (isinstance(dv, float) and math.isnan(dv)):
                continue
            cv = _record_value(rec, k, fid)
            if cv is _UNK:
                continue
            score += 1 if _same(cv, dv) else -2
        if rec.replayed:
            score += 2
        if columns and set(rec.keys) == columns:
            score += 100
        if score > best_score or (score == best_score and best is not None
                                  and (rec.script, rec.line) < (best.script, best.line)):
            best, best_score = rec, score
    return best if best_score >= 3 else None


def _const_kind(rec: CodeRecord, key: str, fid: str) -> tuple[str, str]:
    """('const'|'loopvar'|'expr'|'carry'|'none', code evidence) for one key of one record."""
    node = rec.nodes.get(key)
    if node is None:
        return "none", ""
    where = f"{rec.script}:L{node.lineno}"
    if key in rec.carry_origin:
        return "carry", f"code: {where} `{rec.segs[key]}` copies a column of {rec.carry_origin[key]}"
    if _is_const(node):
        return "const", f"code: {where} literal `{rec.segs[key]}` in the dict built at L{rec.line}"
    if _is_const_ifexp(node):
        env = dict(rec.env)
        env["fid"] = fid
        leaf = _ifexp_leaf(node, env)
        txt = rec.segs[key]
        return "const", f"code: {where} `{txt}` -> branch `{ast.unparse(leaf) if leaf is not None else '?'}`"
    if isinstance(node, ast.Name) and node.id in rec.loop_vars:
        return "loopvar", f"code: {where} `{key}: {rec.segs[key]}` taken from a per-fuel literal tuple"
    loop_deps = [t for t in rec.deps.get(key, frozenset()) if t.split("@")[0] in rec.loop_vars]
    if not loop_deps:
        return "const", f"code: {where} `{rec.segs[key]}` evaluates to the same literal for every fuel of the loop"
    return "expr", f"code: {where} `{rec.segs[key]}`"


def _group_size(ctx: _Ctx, rec: CodeRecord, key: str, value: Any) -> int:
    """How many fuels of this table, built by the same dict literal, carry the same value for key."""
    if ctx.groups is None:
        ctx.groups = {}
        for _, r in ctx.inv.iterrows():
            m = ctx.inv_match.get(str(r["fuel_id"]))
            if m is None:
                continue
            for k in NUMERIC_FIELDS:
                if k in r.index and _num(r[k]) is not None:
                    g = (m.template_id, k, _fmt(r[k]))
                    ctx.groups[g] = ctx.groups.get(g, 0) + 1
    return ctx.groups.get((rec.template_id, key, _fmt(value)), 0)


def _classify_inventory_value(ctx: _Ctx, fid: str, key: str, value: Any) -> tuple[str, str]:
    rec = ctx.inv_match.get(fid)
    kind, code_ev = _const_kind(rec, key, fid) if rec is not None else ("none", "")
    raw = ctx.raw.check_scalar(fid, key, value)
    raw_ev = raw[1] if raw else ""
    rel = ctx.relations.get((fid, key))
    ev = lambda *parts: " | ".join(p for p in parts if p)  # noqa: E731

    if kind == "const":
        if key in STRUCTURAL_ZERO_FIELDS and _close(value, 0.0):
            return "unknown", ev(code_ev, "structural zero literal (no oxygenate / no lump); not a measurement")
        n = _group_size(ctx, rec, key, value)
        if n >= 2:
            note = f"data: the same value {_fmt(value)} is assigned to {n} fuels by this code path"
            if raw and raw[0] == "match":
                note += " (coincides with the raw value for this fuel)"
            return "placeholder_constant", ev(code_ev, note, raw_ev)
        if raw and raw[0] == "match":
            return "raw_file", ev(raw_ev, code_ev)
        if raw and raw[0] == "conflict":
            return "unknown", ev(code_ev, "typed literal conflicts with the raw file", raw_ev)
        return "cited_typed", ev(code_ev, f"single-fuel literal; source label `{_source_of(ctx, fid)}`; "
                                          "not verified against a raw file")
    if kind == "expr" or rel:
        replay_ok = False
        if rec is not None and rec.replayed:
            cv = rec.values.get(key, _UNK)
            replay_ok = cv is not _UNK and _same(cv, value)
        if rel or replay_ok:
            parts = [code_ev, rel or ""]
            if replay_ok:
                parts.append("replaying this expression on the fuel's literal inputs reproduces the stored value")
            if raw and raw[0] == "match":
                parts.append("(also equals the raw value: " + raw_ev + ")")
            elif raw:
                parts.append(raw_ev)
            return "formula_derived", ev(*parts)
    if raw and raw[0] == "match":
        return "raw_file", ev(raw_ev, code_ev)
    if kind == "loopvar":
        if raw and raw[0] == "conflict":
            return "unknown", ev(code_ev, "typed literal conflicts with the raw file", raw_ev)
        return "cited_typed", ev(code_ev, f"source label `{_source_of(ctx, fid)}`; not verified against a raw file",
                                 raw_ev)
    if kind == "carry":
        return "unknown", ev(code_ev, raw_ev, "carried from an intermediate table; see that table's ledger rows")
    return "unknown", ev(code_ev or "no code path reproduces this value", raw_ev)


def _source_of(ctx: _Ctx, fid: str) -> str:
    r = ctx.inv[ctx.inv["fuel_id"] == fid]
    return str(r["source"].iloc[0]) if len(r) and "source" in r.columns else ""


def _classify_claims(ctx: _Ctx, fid: str, claims: str) -> list[tuple[str, str, str, str]]:
    out = []
    for prop in sorted({c.strip() for c in str(claims).split(";") if c.strip()}):
        got = ctx.raw.scalar(fid, prop)
        if got is not None and got[1]:
            path, vals = got
            notes = "; ".join(n for _, n in vals[:3])
            out.append((f"claim:{prop}", "listed", "raw_file",
                        f"claim in available_properties verified: {path} has {prop} for {fid} -> {notes}"))
        elif got is not None:
            out.append((f"claim:{prop}", "listed", "unknown", f"claim NOT found in {got[0]} for {fid}"))
        else:
            out.append((f"claim:{prop}", "listed", "unknown",
                        "claim in available_properties; no raw file maps this fuel/property"))
    return out


# --------------------------------------------------------------------------- compositions
def _comp_records_for(ctx: _Ctx, fid: str, inv_rec: CodeRecord | None) -> list[CodeRecord]:
    """Composition dicts that produced system fuel:fid (replayed for this fuel, else templates of its loop)."""
    sys_id = f"fuel:{fid}"
    scripts = [inv_rec.script] if inv_rec is not None else []
    scripts += sorted({r.script for r in ctx.recs} - set(scripts))
    for script in scripts:
        rep = [r for r in ctx.recs if r.kind == "composition" and r.script == script
               and r.values.get("system_id") == sys_id]
        if rep:
            return rep
    if inv_rec is not None:
        return [r for r in ctx.recs if r.kind == "composition" and r.script == inv_rec.script
                and r.outer_loop_line == inv_rec.outer_loop_line and not r.replayed]
    return []


def _upstream(ctx: _Ctx, fid: str, exclude_script: str) -> tuple[CodeRecord | None, list[CodeRecord]]:
    """For rows a script copied from an earlier table: the upstream script's inventory dict that
    best matches this fuel (by the values both share) and that dict's composition records."""
    inv_row = ctx.inv[ctx.inv["fuel_id"] == fid]
    if not len(inv_row):
        return None, []
    others = [r for r in ctx.recs if r.script != exclude_script]
    up = _match_inventory(fid, inv_row.iloc[0], set(), others)
    if up is None:
        return None, []
    sys_id = f"fuel:{fid}"
    comps = [r for r in others if r.kind == "composition" and r.script == up.script
             and r.values.get("system_id") == sys_id]
    if not comps:
        comps = [r for r in others if r.kind == "composition" and r.script == up.script and not r.replayed
                 and r.outer_loop_line == up.outer_loop_line and (up.outer_loop_line or r.loop_line)]
    return up, comps


def _classify_compositions(ctx: _Ctx) -> list[dict]:
    comp = ctx.comp.copy()
    comp["_fid"] = comp["system_id"].astype(str).str.replace(r"^fuel:", "", regex=True)
    comp["_ord"] = comp.groupby(["system_id", "inchikey"]).cumcount()
    rows: list[dict] = []

    # data relation: identical fraction vectors across systems
    sig = {}
    for sid, g in comp.groupby("system_id", sort=True):
        sig[sid] = tuple(sorted((str(a), round(float(b), 5)) for a, b in zip(g["inchikey"], g["fraction"])))
    groups: dict[tuple, list[str]] = {}
    for sid, s in sig.items():
        groups.setdefault(s, []).append(sid)
    identical = {sid: members for members in groups.values() if len(members) >= 2 for sid in members}

    inv_by_id = {str(r["fuel_id"]): r for _, r in ctx.inv.iterrows()}
    per_system: dict[str, list[tuple[int, str, str]]] = {}

    for sid, g in comp.groupby("system_id", sort=True):
        fid = str(g["_fid"].iloc[0])
        inv_rec = ctx.inv_match.get(fid)
        crecs = _comp_records_for(ctx, fid, inv_rec)
        rep_map: dict[tuple[str, int], CodeRecord] = {}
        counts: dict[str, int] = {}
        for r in crecs:
            if r.replayed:
                # key by InChIKey and by SMILES: v2 looked InChIKeys up in compounds.csv at run time
                # (``smi2ik.get``), which the replay cannot see, so the SMILES is the stable identity.
                for key in (r.values.get("inchikey"), "smi:" + str(r.values.get("smiles", ""))):
                    o = counts.get(key, 0)
                    counts[key] = o + 1
                    rep_map.setdefault((key, o), r)
        g = g.assign(_sord=g.groupby(g["smiles"].fillna("").astype(str)).cumcount())
        raw_comp = ctx.raw.composition(fid)
        gc = ctx.raw.gc_classes(fid)
        results = []
        for idx, row in g.iterrows():
            ik, o, val = str(row["inchikey"]), int(row["_ord"]), row["fraction"]
            rec = rep_map.get((ik, o))
            if rec is None or not _same(rec.values.get("fraction"), val):
                smi = "" if pd.isna(row["smiles"]) else str(row["smiles"])
                alt = rep_map.get(("smi:" + smi, int(row["_sord"]))) if smi else None
                rec = alt if alt is not None else rec
            cls, evid = _classify_one_fraction(ctx, fid, row, ik, o, val, g, rec, crecs,
                                               inv_rec, raw_comp, gc, inv_by_id.get(fid))
            if sid in identical and cls not in ("raw_file", "target_derived_composition"):
                others = [m for m in identical[sid] if m != sid]
                cls = "placeholder_constant"
                evid = (f"data: identical fraction vector across {len(identical[sid])} systems "
                        f"({', '.join(sorted(identical[sid]))}) | " + evid)
            results.append((idx, cls, evid))
        # one composition vector is one object: if any resolved component is target-derived, all are
        if any(c == "target_derived_composition" for _, c, _ in results):
            results = [(i, "target_derived_composition", e if c == "target_derived_composition" else
                        "system vector is target-derived (other components computed from this fuel's own "
                        "properties); this row is its normaliser/remainder | " + e) for i, c, e in results]
        per_system[sid] = results

    for sid, results in per_system.items():
        for idx, cls, evid in results:
            row = comp.loc[idx]
            rows.append({"fuel_id": str(row["_fid"]), "field": f"fraction:{row['inchikey']}#{int(row['_ord'])}",
                         "value": _fmt(row["fraction"]), "provenance_class": cls, "evidence": evid})
    return rows


def _classify_one_fraction(ctx, fid, row, ik, o, val, g, rec, crecs, inv_rec, raw_comp, gc, inv_row):
    smi = "" if pd.isna(row.get("smiles")) else str(row["smiles"]).strip()
    lumped = str(row.get("is_lumped")).lower() == "true"
    basis = str(row.get("basis"))

    # (1) raw molecule list for this fuel: identity (SMILES) + renormalised value + basis
    if raw_comp is not None and not lumped:
        path, rbasis, rdf = raw_comp
        hit = rdf[rdf["smiles"] == smi]
        resolved = g[g["is_lumped"].astype(str).str.lower() != "true"]
        rmap = dict(zip(rdf["smiles"], rdf["frac"]))
        ratios = [float(f) / rmap[s] for s, f in zip(resolved["smiles"].astype(str).str.strip(), resolved["fraction"])
                  if s in rmap and rmap[s] and rmap[s] > 0 and _num(f) is not None]
        if len(hit) == 1 and len(ratios) >= 2 and hit["frac"].iloc[0] > 0:
            k = float(np.median(ratios))
            mine = float(val) / float(hit["frac"].iloc[0])
            if _close(mine, k, rel=3e-3, abs_=0.0) or abs(float(val) - k * hit["frac"].iloc[0]) <= 1e-5:
                note = (f"raw: {path} {smi} frac_value={_fmt(hit['frac'].iloc[0])} ({rbasis}); stored value = raw x "
                        f"{_fmt(round(k, 6))}, the common renormalisation factor of the {len(ratios)} matched "
                        f"components")
                if rbasis.lower().startswith(basis.lower()[:4]):
                    return "raw_file", note
                return "unknown", (note + f" | BASIS MISMATCH: raw frac_type is `{rbasis}` but the table "
                                          f"says basis=`{basis}`")
    if raw_comp is not None and lumped:
        path, rbasis, rdf = raw_comp
        resolved_raw = rdf["frac"].sum() / (100.0 if rdf["frac"].sum() > 1.5 else 1.0)
        if _close(1.0 - resolved_raw, val, rel=1e-3, abs_=1e-3):
            return "raw_file", f"raw: {path} listed components sum to {_fmt(resolved_raw)}; lump = 1 - sum"
        return "cited_typed", (f"lump literal {_fmt(val)} typed in the script; {path} components sum to "
                               f"{_fmt(round(resolved_raw, 4))} (lump would be {_fmt(round(1 - resolved_raw, 4))})")

    # (2) GC hydrocarbon-class table mapped onto representative molecules
    if gc is not None and not lumped:
        path, gdf = gc
        pct = pd.to_numeric(gdf["pct"], errors="coerce")
        totals = gdf["class_or_component"].str.endswith("_total").sum()
        return "formula_derived", (
            f"data: {path} gives {len(gdf)} hydrocarbon-CLASS rows for {fid.replace('_', '-')} (sum of pct = "
            f"{_fmt(round(pct.sum(), 2))}, including sub-totals and `*_total` rows, so classes are double-counted); "
            f"each class pct was rescaled by (1 - hardcoded unresolved) / sum and assigned to one representative "
            f"molecule `{smi}` by a substring rule on the class label -- the molecule identity is not in the raw file")

    # (3) replayed code record for this exact row
    if rec is not None and _same(rec.values.get("fraction"), val):
        rep = ctx.reps[rec.script]
        d = rec.deps.get("fraction", frozenset())
        props = rep.property_deps(d)
        where = f"{rec.script}:L{rec.nodes['fraction'].lineno} `{rec.segs['fraction']}`"
        if props:
            chain = "; ".join(rep.explain(d))
            data_rel = _composition_relation_text(fid, ik, val, inv_row, rec)
            return "target_derived_composition", (
                f"code: {where} depends on this fuel's own {', '.join(props)} via {chain} | "
                f"replaying the code on the fuel's literal inputs reproduces {_fmt(val)}" +
                (f" | {data_rel}" if data_rel else ""))
        node = rec.nodes["fraction"]
        loop_deps = [t for t in d if t.split("@")[0] in rec.loop_vars]
        if _is_const(node) or not loop_deps:
            return "cited_typed", (f"code: {where} literal (no dependence on per-fuel properties); "
                                   f"source label `{_source_of(ctx, fid)}`; not verified against a raw file")
        arg = node.args[0] if isinstance(node, ast.Call) and node.args else node
        lit_tuple = isinstance(arg, ast.Name) and arg.id in rec.loop_vars
        computed = [t for t in d if t in rep.defs and t not in rep.literal_defs]   # arithmetic in the chain
        if lit_tuple and not computed:
            return "cited_typed", (f"code: {where} per-component literal from a recipe table; source label "
                                   f"`{_source_of(ctx, fid)}`; not verified against a raw file")
        return "formula_derived", (f"code: {where} computed from literals and per-fuel inputs "
                                   f"({', '.join(sorted({t.split('@')[0] for t in loop_deps}))}); "
                                   f"replay reproduces {_fmt(val)}")

    # (4) template (non-literal loop): model reconstruction or carry-over
    tmpl = [r for r in crecs if not r.replayed]
    origin = None
    if tmpl:
        t = tmpl[0]
        if "fraction" in t.carry_origin:
            origin = t.carry_origin["fraction"]
        if t.model_calls and not (origin or "").startswith("data/real_fuels"):
            return "formula_derived", (
                f"reconstructed: {t.script}:L{t.nodes['fraction'].lineno} `{t.segs['fraction']}` is produced by "
                f"{', '.join(t.model_calls)} from bulk class contents and D86 points (fuels.csv) with an assumed "
                f"saturate split; not a measured molecule list")
    if origin and origin.startswith("data/real_fuels"):
        # carried from the v2 tables: classify with the code that built that table
        up, upstream = _upstream(ctx, fid, tmpl[0].script)
        if up is not None and upstream:
            match = next((r for r in upstream if r.replayed and _same(r.values.get("fraction"), val) and (
                r.values.get("inchikey") == ik or str(r.values.get("smiles")) == smi)), None)
            sub_cls, sub_ev = _classify_one_fraction(
                ctx, fid, row, ik, o, val, g, match, [r for r in upstream if not r.replayed],
                up, None, None, inv_row)
            sub_ev = sub_ev.replace("reconstructed: ", "", 1)
            prefix = "reconstructed: " if sub_cls == "formula_derived" and \
                any(r.model_calls for r in upstream) else ""
            return sub_cls, f"{prefix}copied from {origin} (upstream dict {up.script}:L{up.line}) | {sub_ev}"
    if origin and origin.startswith("data/processed"):
        return "cited_typed", (f"code: copies `fraction` from processed table {origin} (rows carry a literature "
                               f"reference); not matched against a raw file")
    if lumped:
        return "unknown", "lump/remainder row; no code path or raw file reproduces it"
    return "unknown", (f"no code path or raw file reproduces this fraction" +
                       (f" (loop copies from {origin})" if origin else ""))


def _composition_relation_text(fid, ik, val, inv_row, rec) -> str:
    """Explicit data relation (independent of the replay) when the inventory row carries the inputs."""
    if inv_row is None:
        env = rec.env
        get = lambda k: _num(env.get(k, _UNK)) if env.get(k, _UNK) is not _UNK else None  # noqa: E731
        ron, mon, dens, dcn, unres = get("ron"), get("mon"), get("dens"), get("dcn"), get("unres")
    else:
        ron, mon = _row_num(inv_row, "ron"), _row_num(inv_row, "mon")
        dens, dcn = _row_num(inv_row, "density_kg_m3"), _row_num(inv_row, "dcn")
        unres = _row_num(inv_row, "unresolved_fraction")
    src = "inventory row" if inv_row is not None else "literal tuple in the script"
    if ik in _GASOLINE_ORDER and None not in (ron, mon, dens, unres) and re.fullmatch(r"CRC_FACE_[A-J]", fid):
        f = _face_gasoline_fracs(ron, mon, dens, unres)[_GASOLINE_ORDER.index(ik)]
        if _close(f, val, abs_=1e-5):
            return (f"data: with RON={_fmt(ron)}, MON={_fmt(mon)}, dens={_fmt(dens)} from the {src}, "
                    f"aro=clip((dens-680)/200,0.05,0.45), ole=clip((RON-MON)/20,0.01,0.25), "
                    f"iso=max(0.10,RON/100*(1-aro-ole-unres)), normalised to 1-unres -> {_fmt(f)} == stored")
    if ik in _DIESEL_ORDER and None not in (dcn, dens, unres) and fid.startswith("CRC_FACE_Diesel_"):
        f = _face_diesel_fracs(dcn, dens, unres)[_DIESEL_ORDER.index(ik)]
        if _close(f, val, abs_=1e-5):
            return (f"data: with DCN={_fmt(dcn)}, dens={_fmt(dens)} from the {src}, aro=clip((dens-800)/150), "
                    f"npar=clip((DCN-25)/70), normalised -> {_fmt(f)} == stored")
    return ""


def classify_fuel_rows(inventory: pd.DataFrame, compositions: pd.DataFrame, source_code: dict[str, str],
                       raw_roots=DEFAULT_RAW_ROOTS, repo_root: Path | str = _REPO_ROOT) -> pd.DataFrame:
    """Classify every value of one inventory/composition table pair.

    Returns one row per (fuel_id, field) with columns ``LEDGER_COLUMNS``:
    inventory numeric fields (NaN = no value, skipped), one ``claim:<prop>`` row
    per property listed in ``available_properties``, and one
    ``fraction:<inchikey>#<ordinal>`` row per composition row. Sorted and
    deterministic.
    """
    repo_root = Path(repo_root)
    reps = {s: _replay_cached(s, c) for s, c in sorted(source_code.items())}
    recs = [r for s in sorted(reps) for r in reps[s].records]
    ctx = _Ctx(inv=inventory.reset_index(drop=True), comp=compositions.reset_index(drop=True), recs=recs,
               reps=reps, raw=RawIndex(raw_roots, repo_root), relations=_inventory_relations(inventory))
    cols = set(inventory.columns)
    for _, r in ctx.inv.iterrows():
        ctx.inv_match[str(r["fuel_id"])] = _match_inventory(str(r["fuel_id"]), r, cols, recs)

    out: list[dict] = []
    for _, r in ctx.inv.iterrows():
        fid = str(r["fuel_id"])
        for key in NUMERIC_FIELDS:
            if key not in r.index or _num(r[key]) is None:
                continue
            cls, evid = _classify_inventory_value(ctx, fid, key, r[key])
            out.append({"fuel_id": fid, "field": key, "value": _fmt(r[key]), "provenance_class": cls,
                        "evidence": evid})
        if "available_properties" in r.index and isinstance(r["available_properties"], str):
            for f, v, c, e in _classify_claims(ctx, fid, r["available_properties"]):
                out.append({"fuel_id": fid, "field": f, "value": v, "provenance_class": c, "evidence": e})
    out.extend(_classify_compositions(ctx))
    df = pd.DataFrame(out, columns=LEDGER_COLUMNS)
    df = _apply_component_count_relation(df, ctx)
    return df.sort_values(["fuel_id", "field"], kind="mergesort").reset_index(drop=True)


def _apply_component_count_relation(df: pd.DataFrame, ctx: _Ctx) -> pd.DataFrame:
    """number_of_resolved_components vs the actual composition rows of the system (data relation)."""
    comp = ctx.comp
    counts = comp[comp["is_lumped"].astype(str).str.lower() != "true"].groupby("system_id").size()
    m = df["field"] == "number_of_resolved_components"
    for i in df.index[m]:
        fid = df.at[i, "fuel_id"]
        n_actual = int(counts.get(f"fuel:{fid}", 0))
        claimed = _num(df.at[i, "value"])
        if claimed is not None and int(claimed) != n_actual:
            df.at[i, "evidence"] = (f"data: claims {int(claimed)} resolved components but the composition table "
                                    f"has {n_actual} non-lumped rows for this fuel | " + df.at[i, "evidence"])
            if df.at[i, "provenance_class"] in ("cited_typed", "raw_file"):
                df.at[i, "provenance_class"] = "unknown"
        elif claimed is not None and df.at[i, "provenance_class"] == "unknown":
            df.at[i, "provenance_class"] = "formula_derived"
            df.at[i, "evidence"] = (f"data: equals the {n_actual} non-lumped composition rows of this fuel "
                                    f"(a count derived from the composition table) | " + df.at[i, "evidence"])
    return df
