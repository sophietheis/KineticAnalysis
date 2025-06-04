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

from kinetic_analysis.analysis.analysis_track import (single_track_analysis,
                                     validate_equation)

from kinetic_analysis.utils.utils import read_csv_file

from .app_function import (upload_csv,
                           list_csv_files,
                           browse_directory)


def layout():
    return (
        # Analyse track in vivo tab
        # File selection and display
        dbc.Row([
            html.H4("Upload data",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
            dbc.Col([
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file to analyse "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_analyze_vivo',
                                children=dbc.Button('Select csv file',
                                                    className="mr-2",
                                                    style={"width": "150px"}
                                                    ),
                                multiple=False,
                            ),
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
        # # Choose equation for the analysis
        # dbc.Row([html.H4("Choose equation for the analysis",
        #                  style={"text-align": "center",
        #                         "color": "#10D79B"}),
        #          html.Br(),
        #          dcc.Markdown('''
        #                      By default the model used for analyse track is :
        #                      $$ \\frac{T - x}{cT^2}
        #                      H(T - x), $$
        #                      where $$c$$ is the initiation rate, and $$T$$ is the residence time.
        #
        #                      $$T=M/k$$ where $M$ is the RNA size (aa) and $k$ is the
        #                      elongation rate.
        #                     ''',
        #                       mathjax=True),
        #          html.Div([
        #              html.P(["Equation ",
        #                      html.Span(className="fas fa-question-circle",
        #                                id="faq_equation",
        #                                style={"cursor": "pointer",
        #                                       "marginLeft": "5px"})],
        #                     style={"height": "auto",
        #                            "margin-bottom": "auto"}),
        #              dbc.Tooltip(
        #                  "You can change the equation. Symbols used are"
        #                  " x (fluorescence input), t (elongation rate) and "
        #                  "c (initiation rate). If you need other symbols, "
        #                  "please contact the admin.",
        #                  target="faq_equation"),
        #              dbc.Row([
        #                  dbc.Col([
        #                      dbc.Input(id='equation', type='string', value="",
        #                                style={'width': '300px'}),
        #                      dbc.FormFeedback('Valid equation input!',
        #                                       type='valid'),
        #                      dbc.FormFeedback('Invalid equation input!',
        #                                       type='invalid'),
        #                  ], width=3),
        #                  dbc.Col([
        #                      dbc.Button('Valid equation',
        #                                 id='submit-button-equation',
        #                                 className="mr-2",
        #                                 style={"width": "150px"}, ),
        #                  ], width=3),
        #              ]),
        #          ]),
        #          ]),
        # dbc.Row([
        #     html.Br(),
        #     html.Br(),
        # ]),
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
        dbc.Row([
            html.Div([
                dcc.Checklist(
                    id='checkbox_simu',
                    options=[
                        {'label': 'Tracks from simulation', 'value':
                            'checked'}],
                    value=[],  # default is unchecked
                    labelStyle={'display': 'inline-block',
                                'margin-right': '10px'}
                )
            ])
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
            dcc.Download(id="download-csv")
        ])

    )


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-vivo', 'children'),
        Output("table-container", "children"),
        Output('loading_data_vivo', 'children'),
        Input('browse_directory_analyze_vivo', 'contents'),
    )
    def browse_directory_analyze_vivo(contents):
        df, output =  upload_csv(contents, app)
        if df is None:
            return output, None, None

        table = dash_table.DataTable(
                    data=df.to_dict('records'),
                    columns=[{"name": i, "id": i} for i in
                             df.columns],
                    page_size=15,
                    style_table={'width': '800px', 'overflowX': 'auto'},
        )
        return None, table, None


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
        Output('download-csv', 'data'),
        Input('start-analyze-btn-vivo', 'n_clicks'),
        State('col_track', 'value'),
        State('col_time', 'value'),
        State('col_intensity', 'value'),
        State('dt-param-vivo', 'value'),
        State('prot-length-param-vivo', 'value'),
        State('save-results-name-vivo', 'value'),
        State('checkbox_simu', 'value')
        # State('equation', 'value'),
    )
    def start_analyze_tracks(n_clicks, *params):

        if n_clicks:
            if app.data['csv_to_analyse'] is None:
                return "No CSV file uploaded.", None, None

            try:
                # Read csv file
                df = app.data['csv_to_analyse']
                df.rename(columns={params[0]: 'TRACK_ID',
                                   params[1]: 'FRAME',
                                   params[2]: 'MEAN_INTENSITY_CH1',
                                  },
                         inplace=True)
                dt = float(params[3])
                t = dt / 0.1
                prot_length = float(params[4])
                ids_track = np.unique(df["TRACK_ID"])
                check_simu = 'checked' in params[6]
                
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
                                                   force_analysis=True,
                                                   first_dot=True,
                                                   simulation=check_simu,)
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

                output_path = params[5] + ".csv"
                results.to_csv(output_path, index=False)

                return "Analysis completed and saved successfully!", None, dcc.send_file(output_path)
            except Exception as e:
                return f"Error: {str(e)}", None, None
        raise PreventUpdate


