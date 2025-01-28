import os
import time
import threading
import webview

from threading import Thread

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog

from dash import Dash, html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_spinner

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from .kinetic_function import (single_track_analysis)


def browse_directory(n_clicks, col_name, app):
    if n_clicks:
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder_selected = filedialog.askdirectory()
        root.destroy()
        if folder_selected:
            print(folder_selected)
        else:
            print(None)
        app.data[col_name] = folder_selected
        return f"Directory chosen: {app.data[col_name]}"


def list_csv_files(directory, col_name, app):
    """
    List csv files inside the directory.
    List of files is stored in col_name.
    :param directory:
    :type directory:
    :param col_name:
    :type col_name:
    :param app:
    :type app:
    :return:
    :rtype:
    """
    if app.data[col_name]:
        app.data['csv_files'] = [
            {'label': file, 'value': file}
            for file in os.listdir(app.data[col_name]) if
            file.endswith('.csv')
        ]
        return app.data['csv_files']
    return []