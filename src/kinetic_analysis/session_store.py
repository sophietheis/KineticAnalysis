"""
Per-browser-session storage.

Previously the app kept all uploaded data, selected columns, and settings
in a single `app.data` dict shared by every user/tab hitting the server
(see the removed "Thread safety need to be changed" comment in
app_dash.py). That meant one user's uploaded CSV could overwrite or leak
into another user's session.

This module replaces that with a filesystem-backed cache keyed by a
per-session id. The id itself lives in the browser (a dcc.Store with
storage_type='session', see app_dash.py), so each browser tab gets its
own isolated slice of data, and that slice also survives a page refresh
within the same tab.
"""
import uuid

from flask_caching import Cache

cache = Cache()

# One process-wide dict per session id would defeat the purpose if the
# app is ever run with multiple worker processes (e.g. gunicorn -w 4);
# FileSystemCache stores each session's data on disk instead, so every
# worker process sees the same data for a given session id.
_DEFAULT_CONFIG = {
    "CACHE_TYPE": "FileSystemCache",
    "CACHE_DIR": ".cache",
    "CACHE_DEFAULT_TIMEOUT": 60 * 60 * 4,  # 4 hours
}


def init_cache(app, config=None):
    """Call once, right after creating the Dash app."""
    cache.init_app(app.server, config=config or _DEFAULT_CONFIG)


def new_session_id():
    return str(uuid.uuid4())


def _session_store(session_id):
    if not session_id:
        return {}
    return cache.get(session_id) or {}


def get_session_data(session_id, key, default=None):
    return _session_store(session_id).get(key, default)


def has_session_data(session_id, key):
    return key in _session_store(session_id)


def set_session_data(session_id, key, value):
    if not session_id:
        raise ValueError("Cannot store data without a session id.")
    store = _session_store(session_id)
    store[key] = value
    cache.set(session_id, store)