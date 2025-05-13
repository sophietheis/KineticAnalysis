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

from analysis.analysis_track import (single_track_analysis,
                                     validate_equation)

from utils.utils import read_csv_file

from .app_function import (list_csv_files,
                           browse_directory)


def layout():
    return (
        # Analyse track in vivo tab
        # File selection and display
        dbc.Row([
            dbc.Col([
                html.Div([
                    # Choose a directory
                    html.Label("Choose Directory where your file is"),
                    html.Br(),
                    dbc.Button("Select folder",
                               id="browse_directory_analyze_vivo",
                               className="mr-2",
                               style={"width": "150px"}, ),
                    html.Div(id='directory-analyze-output-vivo', style={
                        'margin-top': '10px'}),
                    html.Br(),
                    # Select a file inside the directory
                    dcc.Dropdown(id='file-dropdown-vivo', options=[],
                                 placeholder="Select a file...",
                                 style={"width": "150px"}, ),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dbc.Button('Validate file',
                                       id='select-file-btn-vivo',
                                       className="mr-2",
                                       style={"width": "150px"}, ),
                        ], width="auto"),
                        dbc.Col([
                            dbc.Spinner(
                                children=[
                                    html.Div(id="loading_data_vivo")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-vivo'),
                ]),
            ], width=3),
            dbc.Col([
                html.Label("DataFrame Visualisation"),
                html.Div(id="table-container", children=[]),
            ], width=5)
        ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Choose columns name
        dbc.Row([
            html.H4("Confirm column name for the analysis",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
            dbc.Col([
                html.Div([
                    html.P("Track ID column",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_track',
                              type='text',
                              value="TRACK_ID",
                              style={'width': '200px'}),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("Time column",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_time',
                              type='text',
                              value="FRAME",
                              style={'width': '200px'}),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("Intensity column",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_intensity',
                              type='text',
                              value="MEAN_INTENSITY_CH1",
                              style={'width': '200px'}),
                ]),
            ]),
        ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Choose equation for the analysis
        dbc.Row([html.H4("Choose equation for the analysis",
                         style={"text-align": "center",
                                "color": "#10D79B"}),
                 html.Br(),
                 dcc.Markdown('''
                             By default the model used for analyse track is : 
                             $$ \\frac{T - x}{cT^2} 
                             H(T - x), $$
                             where $$c$$ is the initiation rate, and $$T$$ is the residence time. 
                             
                             $$T=M/k$$ where $M$ is the RNA size (aa) and $k$ is the 
                             elongation rate. 
                            ''',
                              mathjax=True),
                 html.Div([
                     html.P(["Equation ",
                             html.Span(className="fas fa-question-circle",
                                       id="faq_equation",
                                       style={"cursor": "pointer",
                                              "marginLeft": "5px"})],
                            style={"height": "auto",
                                   "margin-bottom": "auto"}),
                     dbc.Tooltip(
                         "You can change the equation. Symbols used are"
                         " x (fluorescence input), t (elongation rate) and "
                         "c (initiation rate). If you need other symbols, "
                         "please contact the admin.",
                         target="faq_equation"),
                     dbc.Row([
                         dbc.Col([
                             dbc.Input(id='equation', type='string', value="",
                                       style={'width': '300px'}),
                             dbc.FormFeedback('Valid equation input!',
                                              type='valid'),
                             dbc.FormFeedback('Invalid equation input!',
                                              type='invalid'),
                         ], width=3),
                         dbc.Col([
                             dbc.Button('Valid equation',
                                        id='submit-button-equation',
                                        className="mr-2",
                                        style={"width": "150px"}, ),
                         ], width=3),
                     ]),
                 ]),
                 ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Choose parameters for the analysis
        dbc.Row([
            html.H4("Confirm parameters for the analysis",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
            dbc.Col([
                html.Div([
                    html.P("dt (sec)", style={"height": "auto",
                                              "margin-bottom": "auto"}),
                    dcc.Input(id='dt-param-vivo', type='number', value=3),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("Protein length (aa)", style={"height": "auto",
                                                         "margin-bottom": "auto"}),
                    dcc.Input(id='prot-length-param-vivo', type='number',
                              value=800),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("File name to save", style={"height": "auto",
                                                       "margin-bottom": "auto"}),
                    dcc.Input(id='save-results-name-vivo', type='text',
                              value='datas_results'),
                ]),
            ]),
            html.Br(),
        ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Generate Button and Spinner Side by Side
        dbc.Row([
            html.Br(),
            html.Br(),
            dbc.Col([
                dbc.Button('Start Analyse Tracks',
                           id='start-analyze-btn-vivo',
                           className="mr-2",
                           style={"width": "300px"}, ),
            ], width="auto"),

            dbc.Col([
                dbc.Spinner(
                    children=[html.Div(id="loading_analysis_vivo")],
                    size="sm",
                    color="primary",
                    type="border",
                    spinner_style={"margin-left": "10px"}
                )
            ], width="auto"),
            html.Div(id='analyze-output-vivo'),
        ])

    )


def register_callbacks(app):
    @app.callback(
        Output('directory-analyze-output-vivo', 'children'),
        Input('browse_directory_analyze_vivo', 'n_clicks'),
    )
    def browse_directory_analyze_vivo(n_clicks):
        return browse_directory(n_clicks, 'directory_analysis_vivo', app)

    @app.callback(
        Output('file-dropdown-vivo', 'options'),
        Input('directory-analyze-output-vivo', 'children')
    )
    def load_csv_file_vivo(directory):
        return list_csv_files(directory, 'directory_analysis_vivo', app)

    @app.callback(
        Output('selected-file-output-vivo', 'children'),
        Output("table-container", "children"),
        Output('loading_data_vivo', 'children'),
        Input('select-file-btn-vivo', 'n_clicks'),
        State('file-dropdown-vivo', 'value')
    )
    def select_file(n_clicks, selected_file):
        if n_clicks and selected_file:
            app.data['selected_file_vivo'] = selected_file
            df = read_csv_file(
                os.path.join(app.data['directory_analysis_vivo'],
                             app.data['selected_file_vivo']),)

            first_10_rows = df.head(10)

            return (f"You selected: {selected_file}",
                    dash_table.DataTable(
                        data=first_10_rows.to_dict('records'),
                        columns=[{"name": i, "id": i} for i in
                                 first_10_rows.columns]
                    ), None)
        raise PreventUpdate

    @app.callback(
        Output('equation', 'valid'),
        Output('equation', 'invalid'),
        Input('submit-button-equation', 'n_clicks'),
        Input('equation', 'value')
    )
    def validate_input(n_clicks, value):
        if n_clicks :
            if not value :
                return False, False
            v_bool, v_message = validate_equation(value)
            if value and v_bool:
                return True, False
            else:
                return False, True
        return False, False

    @app.callback(
        Output('analyze-output-vivo', 'children'),
        Output('loading_analysis_vivo', 'children'),
        Input('start-analyze-btn-vivo', 'n_clicks'),
        State('file-dropdown-vivo', 'value'),
        State('col_track', 'value'),
        State('col_time', 'value'),
        State('col_intensity', 'value'),
        State('dt-param-vivo', 'value'),
        State('prot-length-param-vivo', 'value'),
        State('save-results-name-vivo', 'value'),
        State('equation', 'value'),
    )
    def start_analyze_tracks(n_clicks, filename, *params):
        if n_clicks:
            try:
                # Read csv file
                df = read_csv_file(os.path.join(app.data[
                                                     'directory_analysis_vivo'],
                                                 filename))
                df.rename(columns={params[0]: 'TRACK_ID',
                                   params[1]: 'FRAME',
                                   params[2]: 'MEAN_INTENSITY_CH1',
                                  },
                         inplace=True)
                dt = float(params[3])
                t = dt / 0.1
                prot_length = float(params[4])
                ids_track = np.unique(df["TRACK_ID"])

                first_time = True
                # Analyse all tracks and save it
                for i in ids_track:
                    datas2 = df[(df["TRACK_ID"] == i)]
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
                                                   normalise_auto=True,
                                                   mm=None,
                                                   rtol=1e-1,
                                                   method="linear",
                                                   force_analysis=False,
                                                   first_dot=False,
                                                   simulation=False,)
                    if first_time:
                        results = pd.DataFrame({"elongation_r": elongation_r,
                                                "init_translation_r": translation_init_r,
                                                "dt": dt,
                                                "id": i, },
                                               index=[0])
                        first_time = False

                    else:
                        results = pd.concat([results,
                                             pd.DataFrame(
                                                 {"elongation_r": elongation_r,
                                                  "init_translation_r": translation_init_r,
                                                  "dt": dt,
                                                  "id": i, }, index=[0])
                                             ], ignore_index=True)

                results.to_csv(
                    os.path.join(app.data['directory_analysis_vivo'],
                                 params[5] + ".csv"))

                return "Analysis completed and saved successfully!", None
            except Exception as e:
                return f"Error: {str(e)}", None
        raise PreventUpdate
