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

from generator.generator_track import (generate_one_track,
                                       generate_tracks,
                                       generate_profile)

from .kinetic_function import (single_track_analysis)
from .app_function import (browse_directory)


def layout():
    return (html.Div([
        # Explanation at the begining of the page
        dbc.Row([
            html.P([
                "In this tab, you will be able to generate tracks according "
                "to a set of parameters. By default, tracks are make with a "
                "time step of 0.1sec. ",
                html.Br(),
                "You can see the fluorescence profile according to the "
                "parameters below. ",
                html.Br(),
                "When happy you can generate a csv file with all the "
                "tracks."]),
            html.Br(),
        ]),
        # Choose directory
        dbc.Row([
            html.Label("Choose Directory to Store Output"),
            html.Br(),
            dbc.Button("Select folder", id="select_directory",
                       className="mr-2",
                       style={"width": "150px"}, ),
            html.Div(id='directory-output', style={'margin-top': '10px'}),
            html.Br(),
        ]),
        # Define parameter of the simulation
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P(["Protein length (aa) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param1",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}
                           ),
                    dbc.Tooltip("Length of the protein in amino acid. This "
                                "value will be added to suntag length.",
                                target="faq_param1"),
                    dcc.Input(id='param1', type='number', value=490,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Suntag length (aa) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param2",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Length of the suntag in amino acid. This "
                                "value will be added to protein length.",
                                target="faq_param2"),
                    dcc.Input(id='param2', type='number', value=796,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Number of suntag ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param3",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Number of suntag repetition",
                                target="faq_param3"),
                    dcc.Input(id='param3', type='number', value=32,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Fluorescence one suntag ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param4",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Fluorescence of one suntag, use as "
                                "reference for fluorescence profile.",
                                target="faq_param4"),
                    dcc.Input(id='param4', type='number', value=4,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Translation rate (aa/sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param5",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Translation rate.",
                                target="faq_param5"),
                    dcc.Input(id='param5', type='number', value=24,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Initiation rate (ribosome/sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param6",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={
                               "height":
                                   "auto",
                               "margin-bottom": "auto"}),
                    dbc.Tooltip("Initiation rate.",
                                target="faq_param6"),
                    dcc.Input(id='param6', type='number', value=1,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Retention time (sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param7",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Protein stay at the RNA for some time",
                                target="faq_param7"),
                    dcc.Input(id='param7', type='number', value=0,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Suntag position (begin or end) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param8",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Choose if Suntag is before or after the "
                                "protein.",
                                target="faq_param8"),
                    dcc.Dropdown(id='param8',
                                 options=[
                                     {'label': 'Begin', 'value': 'begin'},
                                     {'label': 'End', 'value': 'end'}],
                                 value='begin',  # Default value
                                 clearable=False,
                                 # Prevents the user from clearing the selection
                                 style={'width': '200px'}
                                 # Adjust width if needed
                                 ),
                ]),
                html.Div([
                    html.P(["Number of tracks ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param9",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Number of tracks to be generated.",
                                target="faq_param9"),
                    dcc.Input(id='param9', type='number', value=100,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P("File name to save", style={"height": "auto",
                                                       "margin-bottom": "auto"}),
                    dcc.Input(id='param10', type='text', value='datas',
                              style={'width': '200px'}),
                ]),
                html.Br(),
                # Generate Button and Spinner Side by Side
                dbc.Row([
                    dbc.Col([
                        dbc.Button('Start Generate Tracks',
                                   id='start-gen-tracks-btn'),
                    ], width="auto"),

                    dbc.Col([
                        dbc.Spinner(
                            children=[html.Div(id="loading_generate")],
                            size="sm", color="primary", type="border",
                            spinner_style={"margin-left": "10px"}
                        )
                    ], width="auto"),
                ], align="center", style={"margin-top": "10px"}),

                html.Div(id='gen-tracks-output'),
            ], width=3),

            dbc.Col([
                dbc.Button('Show Profile', id='show-profile-btn',
                           className="mr-1"),
                dcc.Graph(id='profile-plot'),

            ], width=5)
        ]),
    ]),
    )


## Callbacks
def register_callbacks(app):
    @app.callback(
        Output("directory-output", "children"),
        [Input("select_directory", "n_clicks")],
    )
    def select_directory(n_clicks):
        """
        Select the directory in which to save the file.
        """
        return browse_directory(n_clicks, 'directory_generation', app)

    @app.callback(
        Output('profile-plot', 'figure'),
        Input('show-profile-btn', 'n_clicks'),
        State('param1', 'value'),
        State('param2', 'value'),
        State('param3', 'value'),
        State('param4', 'value'),
        State('param5', 'value'),
        State('param6', 'value'),
        State('param7', 'value'),
        State('param8', 'value'),
    )
    def update_profile_plot(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """
        if n_clicks:
            try:
                # Generate profile
                x, y = generate_profile(float(params[0]),
                                        float(params[1]),
                                        float(params[2]),
                                        float(params[3]),
                                        float(params[4]),
                                        float(params[6]),
                                        params[7])
                # Create the figure
                figure = make_subplots(rows=3,
                                       cols=1,
                                       subplot_titles=(
                                           'One protein fluo profile',
                                           'One track fluo profile',
                                           'Number of translation'))

                # Plot one protein profile
                figure.add_trace(go.Scatter(x=x, y=y,
                                            mode='lines',
                                            name='Profile one prot'),
                                 row=1,
                                 col=1)
                figure.update_xaxes(title_text='Time', row=1, col=1)
                figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

                # Generate one track
                x, y, y_number = generate_one_track(float(params[0]),
                                                    float(params[1]),
                                                    float(params[2]),
                                                    float(params[3]),
                                                    float(params[4]),
                                                    float(params[5]),
                                                    float(params[6]),
                                                    params[7])

                # Plot one track
                figure.add_trace(go.Scatter(x=x, y=y,
                                            mode='lines',
                                            name='Profile track'),
                                 row=2,
                                 col=1)
                figure.update_xaxes(title_text='Time', row=2, col=1)
                figure.update_yaxes(title_text='Fluorescence', row=2, col=1)

                # Plot number of translation
                figure.add_trace(go.Scatter(x=x, y=y_number,
                                            mode='lines',
                                            name='Profile track'),
                                 row=3,
                                 col=1)
                figure.update_xaxes(title_text='Time', row=3, col=1)
                figure.update_yaxes(title_text='Number of translation', row=3,
                                    col=1)

                figure.update_layout(width=1000, height=800, )
                return figure

            except Exception as e:
                print(e)
                return {
                    'data': [],
                    'layout': go.Layout(title='Error', xaxis={'title': 'Time'},
                                        yaxis={'title': 'Fluorescence'})
                }
        raise PreventUpdate

    @app.callback(
        Output('gen-tracks-output', 'children'),
        Output('loading_generate', 'children'),
        Input('start-gen-tracks-btn', 'n_clicks'),
        State('param1', 'value'),
        State('param2', 'value'),
        State('param3', 'value'),
        State('param4', 'value'),
        State('param5', 'value'),
        State('param6', 'value'),
        State('param7', 'value'),
        State('param8', 'value'),
        State('param9', 'value'),
        State('param10', 'value'),
    )
    def start_generate_tracks(n_clicks, *params):
        # Generate all tracks and save it
        if n_clicks:
            try:
                datas = generate_tracks(int(params[8]),
                                        float(params[0]),
                                        float(params[1]),
                                        float(params[2]),
                                        float(params[3]),
                                        float(params[4]),
                                        float(params[5]),
                                        float(params[6]),
                                        params[7],
                                        )
                datas.to_csv(os.path.join(app.data['directory_generation'],
                                          params[9] + ".csv"))

                return "Tracks generated and saved successfully!", None
            except Exception as e:
                return f"Error: {str(e)}", None
        raise PreventUpdate
