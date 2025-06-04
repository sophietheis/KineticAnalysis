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
                                     validate_equation,
                                     fit_function)
from kinetic_analysis.utils.utils import read_csv_file
from .app_function import (upload_csv,
                           list_csv_files,
                           browse_directory)


def layout():
    return (
        # Analyse track in vivo tab
        # File selection and display
        dbc.Row([
            dbc.Col([
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file to analyse "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_analyze_vivo2',
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
                                    html.Div(id="loading_data_vivo2")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-vivo2'),
                ]),
            ], width=3),
            dbc.Col([
                html.Label("DataFrame Visualisation"),
                html.Div(id="table-container2", children=[]),
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
                    dcc.Input(id='col_track2',
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
                    dcc.Input(id='col_time2',
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
                    dcc.Input(id='col_intensity2',
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
        #                      dbc.Input(id='equation2', type='string', value="",
        #                                style={'width': '300px'}),
        #                      dbc.FormFeedback('Valid equation input!',
        #                                       type='valid'),
        #                      dbc.FormFeedback('Invalid equation input!',
        #                                       type='invalid'),
        #                  ], width=3),
        #                  dbc.Col([
        #                      dbc.Button('Valid equation',
        #                                 id='submit-button-equation2',
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
                    dcc.Input(id='dt-param-vivo2', type='number', value=3),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("Protein length (aa)", style={"height": "auto",
                                                         "margin-bottom": "auto"}),
                    dcc.Input(id='prot-length-param-vivo2', type='number',
                              value=800),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("id of the track to analyse",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='id_track2',
                              type='number',
                              value=0),
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
            dbc.Col([
                dbc.Button('Analyse and show one track',
                           id='analyse_show_button2'),
            ], width="auto"),

            dbc.Col([
                dbc.Spinner(
                    children=[html.Div(id="loading_track_plot")],
                    size="sm", color="primary", type="border",
                    spinner_style={"margin-left": "10px"}
                )
            ], width="auto"),
        ], align="center", style={"margin-top": "10px"}),
        html.Div(id='loading_output1'),
        html.Div(id='loading_output2'),
        dbc.Row([
            dcc.Graph(id='plot_results'),
        ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
    )


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-vivo2', 'children'),
        Output("table-container2", "children"),
        Output('loading_data_vivo2', 'children'),
        Input('browse_directory_analyze_vivo2', 'contents'),
    )
    def browse_directory_analyze_vivo2(contents):
        df, output = upload_csv(contents, app)
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
        Output('equation2', 'valid'),
        Output('equation2', 'invalid'),
        Input('submit-button-equation2', 'n_clicks'),
        State('equation2', 'value')
    )
    def validate_input2(n_clicks, value):
        if n_clicks:
            if not value:
                return False, False
            v_bool, v_message, func_ = validate_equation(value)
            if value and v_bool:
                app.data["equation_f"] = func_
                return True, False
            else:
                return False, True
        return False, False


    @app.callback(
        Output('plot_results', 'figure'),
        Output('loading_output1', 'children'),
        Output('loading_output2', 'children'),
        Output('loading_track_plot', 'children'),
        Input('analyse_show_button2', 'n_clicks'),
        State('col_track2', 'value'),
        State('col_time2', 'value'),
        State('col_intensity2', 'value'),
        State('dt-param-vivo2', 'value'),
        State('prot-length-param-vivo2', 'value'),
        State('id_track2', 'value'),
        # State('equation2', 'value'),
    )
    def analyse_display_track(n_clicks, *params):
        figure = make_subplots(rows=2,
                               cols=1,
                               subplot_titles=["track profile",
                                               "autocorrelation",
                                               ])
        if n_clicks:
            try:
                # Read csv file
                df = app.data['csv_to_analyse']

                df.rename(columns={params[0]: 'TRACK_ID',
                                      params[1]: 'FRAME',
                                      params[2]: 'MEAN_INTENSITY_CH1',
                                      },
                             inplace=True)

                if int(params[5]) not in np.unique(df["TRACK_ID"]):
                    return figure, "This ID does not exist.", "", None

                dt = float(params[3])
                prot_length = float(params[4])
                datas2 = df[(df["TRACK_ID"] == int(params[5]))]
                (x,
                 y,
                 x_auto,
                 y_auto,
                 elongation_r,
                 translation_init_r,
                 perr) = single_track_analysis(datas2,
                                               int(params[5]),
                                               delta_t=dt,
                                               protein_size=prot_length,
                                               normalise_intensity=1,
                                               normalise_auto=True,
                                               mm=None,
                                               rtol=1e-1,
                                               method="linear",
                                               force_analysis=True,
                                               first_dot=False,
                                               simulation=False,
                                               func_ = app.data["equation_f"]
                                               )

                # plot track profile
                figure.add_trace(go.Scatter(x=x,
                                            y=y,
                                            mode="lines",
                                            name="Signal"),
                                            row=1,
                                            col=1,
                                            )
                figure.update_xaxes(title_text='Time (sec)', row=1, col=1)
                figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

                # plot track profile
                figure.add_trace(go.Scatter(x=x_auto,
                                            y=y_auto,
                                            mode="lines",
                                            name="Autocorrelation"),
                                            row=2,
                                            col=1),

                figure.add_trace(go.Scatter(x=x_auto[:int(len(x_auto)/2)],
                                            y=fit_function(x_auto,
                                                           prot_length/elongation_r, 1/translation_init_r)[:int(len(x_auto)/2)],
                                            mode="lines",
                                            name="Fit"),
                                 row=2,
                                 col=1),

                figure.update_xaxes(title_text='Time delay (tau)',
                                    row=2,
                                    col=1)
                figure.update_yaxes(title_text='G(tau)', row=2, col=1)
                figure.update_layout(width=1000, height=800, )

                str_output1 = (f"elongation rate: "
                              f"{elongation_r:.2f} aa/sec ")
                str_output2 = (f"initiation rate: "
                              f"{translation_init_r:.2f} rib/sec")

                return figure, str_output1, str_output2, None
            except Exception as e:
                print(e)
                return {
                    'data': [],
                    'layout': go.Layout(title='Error', xaxis={'title': 'Time'},
                                        yaxis={'title': 'Fluorescence'})
                }, str(e),"",  None
            raise PreventUpdate
        return figure, "", "", None
