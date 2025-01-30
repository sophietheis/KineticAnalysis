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

from analysis.analysis_track import (single_track_analysis)

from .app_function import (list_csv_files,
                          browse_directory)

def layout():
    return(
        # Analyse track simulation tab

            html.Div([

                # Choose a directory
                html.Label("Choose Directory where your file is"),
                html.Br(),
                dbc.Button("Select folder", id="browse_directory_analyze",
                           className="mr-2", style={"width": "150px"}, ),
                html.Div(id='directory-analyze-output',
                         style={'margin-top': '10px'}),
                html.Br(),
                # Select a file inside the directory
                dcc.Dropdown(id='file-dropdown', options=[],
                             placeholder="Select a file...",
                             style={"width": "150px"}, ),
                dbc.Button('Validate file', id='select-file-btn',
                           className="mr-2", style={"width": "150px"}, ),
                html.Div(id='selected-file-output'),
                html.Br(),
                html.Br(),
                html.Br(),
                dbc.Col([
                    html.Div([
                        html.P("dt", style={"height": "auto",
                                            "margin-bottom": "auto"}),
                        dcc.Input(id='dt-param', type='number', value=3),
                    ]),
                    html.Div([
                        html.P("Protein length (aa)", style={"height": "auto",
                                                             "margin-bottom": "auto"}),
                        dcc.Input(id='prot-length-param', type='number',
                                  value=800),
                    ]),
                    html.Div([
                        html.P("File name to save", style={"height": "auto",
                                                           "margin-bottom": "auto"}),
                        dcc.Input(id='save-results-name', type='text',
                                  value='datas_results'),
                    ]),
                    html.Br(),
                    # Generate Button and Spinner Side by Side
                    dbc.Row([
                        dbc.Col([
                            dbc.Button('Start Analyze Tracks',
                                       id='start-analyze-btn',
                                       className="mr-2",
                                       style={"width": "150px"}, ),
                        ], width="auto"),

                        dbc.Col([
                            dbc.Spinner(
                                children=[html.Div(id="loading_analysis")],
                                size="sm", color="primary", type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                        html.Div(id='analyze-output'),
                    ], align="center", style={"margin-top": "10px"}),
                ], width=3),
            ]),

    )


def register_callbacks(app):
    @app.callback(
        Output('directory-analyze-output', 'children'),
        Input('browse_directory_analyze', 'n_clicks'),
    )
    def browse_directory_analyze(n_clicks):
        return browse_directory(n_clicks, 'directory_analysis', app)

    @app.callback(
        Output('file-dropdown', 'options'),
        Input('directory-analyze-output', 'children')
    )
    def load_csv_file_simu(directory):
        return list_csv_files(directory, 'directory_analysis', app)

    @app.callback(
        Output('selected-file-output', 'children'),
        Input('select-file-btn', 'n_clicks'),
        State('file-dropdown', 'value')
    )
    def select_file(n_clicks, selected_file):
        if n_clicks and selected_file:
            app.data['selected_file'] = selected_file
            return f"You selected: {selected_file}"
        raise PreventUpdate

    @app.callback(
        Output('analyze-output', 'children'),
        Output('loading_analysis', 'children'),
        Input('start-analyze-btn', 'n_clicks'),
        State('file-dropdown', 'value'),
        State('dt-param', 'value'),
        State('prot-length-param', 'value'),
        State('save-results-name', 'value'),
    )
    def start_analyze_tracks(n_clicks, filename, *params):
        if n_clicks:
            try:
                # Read csv file
                datas = pd.read_csv(os.path.join(app.data['directory_analysis'], filename),
                                    index_col="Unnamed: 0")
                dt = float(params[0])
                t = dt / 0.1
                prot_length = float(params[1])
                nb_track = len(np.unique(datas["TRACK_ID"]))

                first_time = True
                # Analyse all tracks and save it
                for i in range(nb_track):
                    datas2 = datas[(datas["TRACK_ID"] == i)][::int(t)]

                    (x,
                     y,
                     x_auto,
                     y_auto,
                     elongation_r,
                     translation_init_r,
                     perr) = single_track_analysis(datas2,
                                                   i,
                                                   delta_t=dt,
                                                   protein_size=prot_length,
                                                   normalise_intensity=1,
                                                   normalize_auto=True,
                                                   mm=None,
                                                   lowpass_=False,
                                                   cutoff=100,
                                                   rtol=1e-1,
                                                   method="linear",
                                                   force_analysis=True,
                                                   first_dot=True,
                                                   simulation=True)
                    if first_time:
                        results = pd.DataFrame({"elongation_r": elongation_r,
                                                "init_translation_r": translation_init_r,
                                                "dt": dt,
                                                "id": i, },
                                               index=[0])
                        first_time = False

                    else:
                        results = pd.concat([results,
                                             pd.DataFrame({"elongation_r": elongation_r,
                                                           "init_translation_r": translation_init_r,
                                                           "dt": dt,
                                                           "id": i, }, index=[0])
                                             ], ignore_index=True)


                results.to_csv(os.path.join(app.data['directory_analysis'], params[2] + ".csv"))

                return "Analysis completed and saved successfully!", None
            except Exception as e:
                return f"Error: {str(e)}", None
        raise PreventUpdate
