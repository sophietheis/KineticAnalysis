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

from kinetic_analysis.generator.generator_track import (generate_one_track,
                                       generate_tracks,
                                       generate_profile)

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
                                      id="faq_param_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}
                           ),
                    dbc.Tooltip("Length of the protein in amino acid. This "
                                "value will be added to suntag length.",
                                target="faq_param_prot_length"),
                    dcc.Input(id='param_prot_length', type='number', value=490,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Suntag length (aa) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_suntag_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Length of the suntag in amino acid. This "
                                "value will be added to protein length.",
                                target="faq_param_suntag_length"),
                    dcc.Input(id='param_suntag_length', type='number', value=796,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Number of suntag ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_nb_suntag",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Number of suntag repetition",
                                target="faq_param_nb_suntag"),
                    dcc.Input(id='param_nb_suntag', type='number', value=32,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Fluorescence one suntag ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_fluo_one_suntag",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Fluorescence of one suntag, use as "
                                "reference for fluorescence profile.",
                                target="faq_param_fluo_one_suntag"),
                    dcc.Input(id='param_fluo_one_suntag', type='number', value=4,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Translation rate (aa/sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_translation_rate",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Translation rate.",
                                target="faq_param_translation_rate"),
                    dcc.Input(id='param_translation_rate', type='number', value=24,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Initiation rate (ribosome/sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_initiation_rate",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={
                               "height":
                                   "auto",
                               "margin-bottom": "auto"}),
                    dbc.Tooltip("Initiation rate.",
                                target="faq_param_initiation_rate"),
                    dcc.Input(id='param_initiation_rate', type='number', value=1,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Retention time (sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_retention_time",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Protein stay at the RNA for some time",
                                target="faq_param_retention_time"),
                    dcc.Input(id='param_retention_time', type='number', value=0,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Suntag position (begin or end) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_pos_suntag",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Choose if Suntag is before or after the "
                                "protein.",
                                target="faq_param_pos_suntag"),
                    dcc.Dropdown(id='param_pos_suntag',
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
                    html.P(["Noise",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_noise",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Add random noise to the protein "
                                "fluorescence profile. 0 if no noise. ",
                                target="faq_param_noise"),
                    dcc.Input(id='param_noise', type='number', value=0.,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["dt (sec)",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_dt",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(
                        "Time step use to generate track. This value "
                        "should be quite small regarding the time step "
                        "used for the analysis",
                        target="faq_param_dt"),
                    dcc.Input(id='param_dt', type='number', value=0.1,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Length of one track (sec)",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(
                        "Length duration of one track."
                        "used for the analysis",
                        target="faq_param_length"),
                    dcc.Input(id='param_length', type='number', value=6000,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Number of tracks ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_nb_tracks",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Number of tracks to be generated.",
                                target="faq_param_nb_tracks"),
                    dcc.Input(id='param_nb_tracks', type='number', value=100,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P("File name to save", style={"height": "auto",
                                                       "margin-bottom": "auto"}),
                    dcc.Input(id='param_filename', type='text', value='datas',
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

            # Show plot profile
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
        State('param_prot_length', 'value'),
        State('param_suntag_length', 'value'),
        State('param_nb_suntag', 'value'),
        State('param_fluo_one_suntag', 'value'),
        State('param_translation_rate', 'value'),
        State('param_initiation_rate', 'value'),
        State('param_retention_time', 'value'),
        State('param_pos_suntag', 'value'),
        State('param_noise', 'value'),
        State('param_dt', 'value'),
        State('param_length', 'value'),
    )
    def update_profile_plot(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """
        if n_clicks:
            try:
                # Generate profile
                noise = False
                if params[8] > 0:
                    noise = True
                x, y = generate_profile(prot_length=float(params[0]),
                                        suntag_length=float(params[1]),
                                        nb_suntag=float(params[2]),
                                        fluo_one_suntag=float(params[3]),
                                        translation_rate=float(params[4]),
                                        retention_time=float(params[6]),
                                        suntag_pos=params[7],
                                        noise=noise,
                                        noise_std=float(params[8]),
                                        step=float(params[9]))
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
                figure.update_xaxes(title_text='Time (sec)', row=1, col=1)
                figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

                # Generate one track
                x, y, y_number = generate_one_track(prot_length=float(params[0]),
                                                    suntag_length=float(params[1]),
                                                    nb_suntag=float(params[2]),
                                                    fluo_one_suntag=float(params[3]),
                                                    translation_rate=float(params[4]),
                                                    binding_rate=float(params[5]),
                                                    retention_time=float(params[6]),
                                                    suntag_pos=params[7],
                                                    noise=noise,
                                                    noise_std=float(params[8]),
                                                    step=float(params[9]),
                                                    length=float(params[10])
                                                    )

                # Plot one track
                figure.add_trace(go.Scatter(x=x, y=y,
                                            mode='lines',
                                            name='Profile track'),
                                 row=2,
                                 col=1)
                figure.update_xaxes(title_text='Time (sec)', row=2, col=1)
                figure.update_yaxes(title_text='Fluorescence', row=2, col=1)

                # Plot number of translation
                figure.add_trace(go.Scatter(x=x, y=y_number,
                                            mode='lines',
                                            name='Profile track'),
                                 row=3,
                                 col=1)
                figure.update_xaxes(title_text='Time (sec)', row=3, col=1)
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
        State('param_prot_length', 'value'),
        State('param_suntag_length', 'value'),
        State('param_nb_suntag', 'value'),
        State('param_fluo_one_suntag', 'value'),
        State('param_translation_rate', 'value'),
        State('param_initiation_rate', 'value'),
        State('param_retention_time', 'value'),
        State('param_pos_suntag', 'value'),
        State('param_noise', 'value'),
        State('param_dt', 'value'),
        State('param_length', 'value'),
        State('param_nb_tracks', 'value'),
        State('param_filename', 'value'),
    )
    def start_generate_tracks(n_clicks, *params):
        # Generate all tracks and save it
        if n_clicks:
            try:
                noise = False
                if params[8] > 0:
                    noise = True
                datas = generate_tracks(n=int(params[11]),
                                        prot_length=float(params[0]),
                                        suntag_length=float(params[1]),
                                        nb_suntag=float(params[2]),
                                        fluo_one_suntag=float(params[3]),
                                        translation_rate=float(params[4]),
                                        binding_rate=float(params[5]),
                                        retention_time=float(params[6]),
                                        suntag_pos=params[7],
                                        noise=noise,
                                        noise_std=float(params[8]),
                                        step=float(params[9]),
                                        length=float(params[10]),
                                        )
                datas.to_csv(os.path.join(app.data['directory_generation'],
                                          params[12] + ".csv"))

                return "Tracks generated and saved successfully!", None
            except Exception as e:
                return f"Error: {str(e)}", None
        raise PreventUpdate
