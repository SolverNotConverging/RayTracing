"""Preserve independent analytical transport checks removed during cleanup."""
import importlib.util
import json
from pathlib import Path


def test_beam_analytical_reference_suite(tmp_path):
    path=Path(__file__).resolve().parents[2]/"benchmarks/beam_transport.py"
    spec=importlib.util.spec_from_file_location("analytical_beam_validation",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.run(tmp_path)
    summary=json.loads((tmp_path/"summary.json").read_text())
    assert summary["passed"]
