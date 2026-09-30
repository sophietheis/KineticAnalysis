import logging


def configure_logging(level=logging.INFO):
    """
    Configure root logging once for the whole app.

    Call this from the entry point (app_dash.py) only. Every other module
    should just do `logger = logging.getLogger(__name__)` and log through
    that, instead of using print().
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )