"""HANNA (Specht et al. 2024, MIT licence) binary activity coefficients, pretrained, no refit.

Needs the cloned repo (code, `Model/HANNA_Val.pt`, `Model/scalerHANNA_Val.pkl`) and the ChemBERTa-77M-MTR weights.
Locations: env `HANNA_REPO` (default scratch/hanna_dl/repo) and `HF_HOME` (default scratch/hanna_dl/hf). Binary mixtures only.
The same call shape as the UNIFAC gamma function in `baselines/flash_unifac.py`: `make_hanna_gamma_fn(smiles)` returns `f(T_K, x) -> gamma`.
"""
from __future__ import annotations

import importlib
import os
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
_STATE: dict = {}


def _paths():
    repo = Path(os.environ.get("HANNA_REPO", ROOT / "scratch/hanna_dl/repo"))
    os.environ.setdefault("HF_HOME", str(ROOT / "scratch/hanna_dl/hf"))
    return repo


def _load():
    if _STATE:
        return _STATE
    import torch
    from rdkit import Chem  # noqa: F401
    repo = _paths()
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    H = importlib.import_module("utils.HANNA")
    U = importlib.import_module("utils.Utils")
    importlib.import_module("utils.Own_Scaler")          # the scaler pickle refers to utils.Own_Scaler.CustomScaler
    model = H.HANNA()
    model.load_state_dict(torch.load(repo / "Model/HANNA_Val.pt", map_location="cpu"))
    model.eval()
    with open(repo / "Model/scalerHANNA_Val.pkl", "rb") as f:
        scaler = pickle.load(f)
    from transformers import AutoModel
    from tokenizers import Regex, Tokenizer
    from tokenizers.models import WordLevel
    from tokenizers.pre_tokenizers import Split
    name = "DeepChem/ChemBERTa-77M-MTR"
    from huggingface_hub import hf_hub_download
    vocab = hf_hub_download(name, "vocab.json")          # transformers 5 no longer writes vocab.json via save_vocabulary
    ctok = Tokenizer(WordLevel.from_file(vocab, unk_token="[UNK]"))
    ctok.pre_tokenizer = Split(pattern=Regex(r"\[(.*?)\]|Br|Cl|."), behavior="isolated")
    _STATE.update(model=model, scaler=scaler, U=U, bert=AutoModel.from_pretrained(name).eval(), ctok=ctok, emb={})
    return _STATE


def _emb(smi: str) -> np.ndarray:
    s = _load()
    can = s["U"].canonicalize_smiles(smi)
    if can not in s["emb"]:
        s["emb"][can] = s["U"].get_smiles_embedding(can, s["ctok"], s["bert"], "cpu").flatten()
    return s["emb"][can]


def ln_gamma_hanna(smiles, x, T_K: float) -> np.ndarray:
    """ln gamma of both components of a binary mixture at mole fractions `x` (x1 clipped to [0, 1], renormalised)."""
    if len(smiles) != 2:
        raise NotImplementedError("HANNA is binary-only")
    import torch
    s = _load()
    x = np.asarray(x, float)
    x1 = float(x[0] / x.sum())
    row = np.concatenate(([T_K, x1], _emb(smiles[0]), _emb(smiles[1])))[None, :]
    X = s["U"].preprocess_input(row)
    T, xx, FP = s["U"].split_and_reshape_input(s["scaler"].transform(X))
    ln_g, _ = s["model"](torch.tensor(T, dtype=torch.float32), torch.tensor(xx, dtype=torch.float32), torch.tensor(FP, dtype=torch.float32))
    return ln_g.detach().numpy()[0]


def gamma_hanna(smiles, x, T_K: float) -> np.ndarray:
    return np.exp(ln_gamma_hanna(smiles, x, T_K))


def make_hanna_gamma_fn(smiles):
    """Gamma function `f(T_K, x)` in the shape used by `predict_flash_le_chatelier`; raises NotImplementedError for more than 2 components."""
    if len(smiles) != 2:
        raise NotImplementedError("HANNA is binary-only")
    return lambda T_K, x: gamma_hanna(smiles, x, float(T_K))
