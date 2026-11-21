import pytest

torch = pytest.importorskip("torch")
from fuelblend.ml import dmpnn  # noqa: E402


def test_device_resolution_reports_actual(monkeypatch):
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    monkeypatch.setenv("FUELBLEND_DEVICE", "cuda")
    with pytest.raises(RuntimeError, match="CUDA requested"):
        dmpnn.resolve_device()
    monkeypatch.delenv("FUELBLEND_DEVICE")
    assert dmpnn.resolve_device().type == "cpu"
    assert dmpnn.resolve_device("cpu").type == "cpu"


def test_device_info_keys(monkeypatch):
    monkeypatch.delenv("FUELBLEND_DEVICE", raising=False)
    info = dmpnn.device_info()
    assert {"device", "torch_version", "cuda_version", "gpu_name"} <= set(info)
    assert info["device"] == "cpu"
