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

from tab_generate_track import layout as tab1_layout
from tab_generate_track import register_callbacks as tab1_callbacks

from tab_analyse_simu import layout as tab2_layout
from tab_analyse_simu import register_callbacks as tab2_callbacks

from tab_analyse_invivo import layout as tab3_layout
from tab_analyse_invivo import register_callbacks as tab3_callbacks


app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
app.title = "Kinetic analysis app"

# Global variables to store states
app.data = {
    'directory_generation': None,
    'directory_analysis': None,
    'directory_analysis_vivo': None,
    'csv_files': [],
    'fig': None,
    'selected_file': None,
}


app.layout = dbc.Container([
    html.H1("Kinetic Analysis"),
    html.P("This tool is made to estimate the kinetic parameters of translation"
           " (initiation rate and elongation rate). "),
    html.P("You can find different tabs: "),
    html.Li("Track generator: generate tracks to test parameters"),
    html.Li("Track analysis simulation : analyse tracks from simulation"),
    html.Li("Track analysis in vivo : analyse tracks from in vivo "
            "experiments"),
    html.Br(),

    dbc.Tabs(id='tabs', active_tab='tab-1', children=[
        dbc.Tab(label="Generate tracks", tab_id="tab-1"),
        dbc.Tab(label="Analyse tracks simu", tab_id="tab-2"),
        dbc.Tab(label="Analyse tracks in vivo", tab_id="tab-3"),
    ]),
    html.Div(id='tabs-content')
])

@app.callback(
    Output('tabs-content', 'children'),
    Input('tabs', 'active_tab')
)
def render_content(tab):
    if tab == 'tab-1':
        return tab1_layout()
    elif tab == 'tab-2':
        return tab2_layout()
    elif tab == 'tab-3':
        return tab3_layout()

# Register callbacks
tab1_callbacks(app)
tab2_callbacks(app)
tab3_callbacks(app)




def run_app():
    app.run_server(debug=False, port=8080)

if __name__ == '__main__':
    # app.run_server(debug=True, port=8080)
    t = Thread(target=run_app)
    t.daemon = True
    t.start()

    window = webview.create_window('Kinetic analysis', 'http://127.0.0.1:8080/',
                                   width=1800, height=1000,
)
    webview.start(debug=False)

