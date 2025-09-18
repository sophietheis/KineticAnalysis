# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['app_dash.py'],
    pathex=[],
    binaries=[],
    datas=[("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/dash", "dash"),
	    ("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/dash_html_components", "dash_html_components"),
	    ("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/dash_core_components", "dash_core_components"),
	    ("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/dash_bootstrap_components", "dash_bootstrap_components"),
	    ("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/dash_spinner", "dash_spinner"),
	    ("/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/dash_table", "dash_table"),],
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
