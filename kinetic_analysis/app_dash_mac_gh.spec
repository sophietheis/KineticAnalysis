# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['app_dash.py'],
    pathex=[],
    binaries=[],
    datas=[("c:/hostedtoolcache/windows/python/3.9.13/x64/lib/site-packages/dash", "dash"),
	    ("c:/hostedtoolcache/windows/python/3.9.13/x64/lib/site-packages/dash_html_components", "dash_html_components"),
	    ("c:/hostedtoolcache/windows/python/3.9.13/x64/lib/site-packages/dash_core_components", "dash_core_components"),
	    ("c:/hostedtoolcache/windows/python/3.9.13/x64/lib/site-packages/dash_bootstrap_components", "dash_bootstrap_components"),
	    ("c:/hostedtoolcache/windows/python/3.9.13/x64/lib/site-packages/dash_spinner", "dash_spinner"),
	    ("c:/hostedtoolcache/windows/python/3.9.13/x64/lib/site-packages/dash_table", "dash_table"),],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='KineticAnalysisMac',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
