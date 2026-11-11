"""Load the project configuration (``config.yaml`` at the repository root).

This is the one place that knows where the config file lives, so every
script resolves data paths the same way.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = REPO_ROOT / "config.yaml"


def load_config(path: str | Path = DEFAULT_CONFIG_PATH) -> dict[str, Any]:
    """Read the YAML config and return it as a nested dict."""
    with open(path) as fh:
        return yaml.safe_load(fh)


def project_path(config: dict[str, Any], key: str) -> Path:
    """Return the absolute path for ``config["paths"][key]``.

    Paths in the config are relative to the repository root.
    """
    return REPO_ROOT / config["paths"][key]
