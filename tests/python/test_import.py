"""Exercise IDE source-root imports in fresh interpreters."""
import os
from pathlib import Path
import subprocess
import sys


SOURCE = Path(__file__).resolve().parents[2] / "python"


def test_source_root_uses_installed_native_extension():
    env = {**os.environ, "PYTHONPATH": str(SOURCE)}
    result = subprocess.run(
        [sys.executable, "-c", "import raytracing as rt; "
         "from pathlib import Path; "
         "assert Path(rt.__file__).parent == Path(" + repr(str(SOURCE / "raytracing")) + "); "
         "assert rt.example_pattern().is_file(); rt.Scene()"],
        env=env, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_missing_native_extension_has_actionable_error():
    env = {**os.environ, "PYTHONPATH": str(SOURCE)}
    result = subprocess.run(
        [sys.executable, "-S", "-c", "import raytracing"],
        env=env, capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "native _core extension is missing" in result.stderr
    assert "partially initialized" not in result.stderr
