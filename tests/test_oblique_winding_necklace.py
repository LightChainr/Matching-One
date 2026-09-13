"""Expose the focused geometric controls to unittest discovery."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from oblique_winding_necklace import GeometryTests  # noqa: E402,F401
