"""
texts.py
========
Loads all user-facing text from a JSON locale file.

To switch language, set the APP_LANG environment variable (default: "en").
The JSON files live in the locales/ folder, next to the utils/ folder.

Usage
-----
    from kineticanalysis.utils.texts import t

    dcc.Markdown(t("equations.msd"), mathjax=True)
    dbc.Tooltip(t("tooltips.prot_length"), target="faq_param_prot_length")
    html.H4(t("labels.section_upload"))
    html.P(t("tab_intros.msd"))

    # Strings with {placeholders} (do not use on LaTeX strings)
    t("misc.some_key", n=3)
"""

import json
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Language selection
# Set the APP_LANG environment variable, e.g. APP_LANG=fr -> locales/fr.json
# ---------------------------------------------------------------------------
LANGUAGE = os.environ.get("APP_LANG", "en")

# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------
_LOCALES_DIR = Path(__file__).parent.parent / "locales"


def _load(language: str) -> dict:
    path = _LOCALES_DIR / f"{language}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"Locale file not found: {path}\n"
            f"Available locales: {[f.stem for f in _LOCALES_DIR.glob('*.json')]}"
        )
    with open(path, encoding="utf-8") as f:
        return json.load(f)


T = _load(LANGUAGE)


def t(key: str, **fmt) -> str:
    """Return the text for a dotted key, e.g. t("labels.prot_length").

    If keyword arguments are given, they are passed to str.format().
    Raises a KeyError naming the full key if it does not exist.
    """
    node = T
    for part in key.split("."):
        try:
            node = node[part]
        except (KeyError, TypeError):
            raise KeyError(f"Missing text key: '{key}'") from None
    return node.format(**fmt) if fmt else node


# Convenience aliases, so existing code using them keeps working.
EQUATIONS = T["equations"]
TOOLTIPS = T["tooltips"]
LABELS = T["labels"]
TAB_INTROS = T["tab_intros"]
MISC = T["misc"]
