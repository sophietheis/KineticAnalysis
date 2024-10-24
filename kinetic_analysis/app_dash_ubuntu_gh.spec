# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['app_dash.py'],
    pathex=[],
    binaries=[],
    datas=[("/opt/hostedtoolcache/Python/3.11.10/x64/lib/python3.11/site-packages/dash", "dash"),
	    ("/opt/hostedtoolcache/Python/3.11.10/x64/lib/python3.11/site-packages/dash_html_components", "dash_html_components"),
	    ("/opt/hostedtoolcache/Python/3.11.10/x64/lib/python3.11/site-packages/dash_core_components", "dash_core_components"),
	    ("/opt/hostedtoolcache/Python/3.11.10/x64/lib/python3.11/site-packages/dash_bootstrap_components", "dash_bootstrap_components"),
	    ("/opt/hostedtoolcache/Python/3.11.10/x64/lib/python3.11/site-packages/dash_spinner", "dash_spinner"),
	    ("/opt/hostedtoolcache/Python/3.11.10/x64/lib/python3.11/site-packages/dash_table", "dash_table"),],
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
    name='KineticAnalysisUbuntu',
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
