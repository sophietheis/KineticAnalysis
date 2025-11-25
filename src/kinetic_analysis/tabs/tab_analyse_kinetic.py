import numpy as np
import pandas as pd

from dash import  html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from .app_function import (upload_csv)
from ..analysis.analysis_track import (single_track_analysis,
                                       validate_equation,
                                       fit_function,
                                       fit_function_exact,
                                       check_track_validity)


def layout():
    items_equation = ["",
                      "((t - x) / (c * t ** 2)) * heaviside((t - x), 0)",
                      ]
    items_solver = ["",
                    "Exact equation",
                    "Approximate equation",
                    "Approximate epitope",
                    "linear fit"]

    return (
        # File selection and display
        html.H4(children="Upload data",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),
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
                                size="lm", #"sm"
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
            ], width=9)
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
                    html.P(["Time column",
                            html.Span(className="fas fa-question-circle",
                                       id="faq_time_col",
                                       style={"cursor": "pointer",
                                              "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_time2',
                              type='text',
                              value="FRAME",
                              style={'width': '200px'}),
                    dbc.Tooltip(
                        "It corresponds to a FRAME column. This will be "
                        "multiply by the dt time step.",
                        target="faq_time_col"),
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
        # Choose equation for the analysis
        dbc.Row([
            html.H4("Choose equation and solver for the analysis",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
            # Solver
            dbc.Col([
                html.H5("Solver"),
                dcc.Markdown('''
                     Choose the solver. \n
                     You can choose between exact_fit, approximation_fit 
                     and linear_fit. \n 
                    ''',
                             mathjax=True),

                dbc.Select(options=items_solver,
                           id="choose-solver1",
                           className="mr-2",
                           style={"width":350},),

                dcc.Markdown('''
                            If you choose exact_fit, the default the 
                             equation used is :
                             $$ G(T) = \\frac{k}{c}(\\frac{2}{3})\\frac{
                             1}{(N(N+1))^2}e^{-kT} \\sum_{n=0}^N(N-n)(
                             N-n+1)(2N+n+1)\\frac{(kT)^n}{n!} $$
                             where $$c$$ is the initiation rate, $$k$$ 
                             the elongation rate and $$T$$ is 
                             the residence time.

                             If you choose approximation_fit, the default the 
                             equation used is :
                             $$ G(x) = \\frac{T - x}{cT^2}
                             H(T - x), $$
                             where $$c$$ is the initiation rate, and $$T$$ is the residence time.

                             $$T=M/k$$ where $M$ is the RNA size (aa) and $k$ is the
                             elongation rate.
                            ''',
                             mathjax=True),
            ]),

            # Equation
            dbc.Col([
                html.H5("Equation"),

                dcc.Markdown('''
                You can change the equation. Symbols used are: 
                * x (fluorescence input), 
                * t (elongation rate) 
                * c (initiation rate). 
                
                If you need other symbols, please contact the admin.
                '''),
                dcc.Markdown(''' This part is not working (and it is not 
                useful) ''',
                             style={'color': 'red', }),
                dbc.Select(options=items_equation,
                            id="equation1",
                            className="mr-2",
                            style={"width": 350}, ),

                html.Br(),


                dbc.Input(id='equation2',
                           type='string',
                           value="",
                           placeholder="Enter equation",
                           style={'width': '350px'}),
                dbc.FormFeedback('Valid equation input!',
                                  type='valid'),
                dbc.FormFeedback('Invalid equation input!',
                                  type='invalid'),

                html.Br(),
                dbc.Button('Validate equation',
                            id='submit-button-equation2',
                            className="mr-2",
                            style={"width": "150px"}, ),


                # display the chosen equation and that it is valid
                dbc.Row([
                     html.Div(id='loading_equation_output2'),
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
        ]),
        # html.Br(),
        # Put here the chosen equation
        dbc.Row([
            dbc.Col(html.P("The chosen solver is: ", className="mb-0"),
                    width="auto"),
            dbc.Col(html.Div(id="choosen-solver1"), width="auto"),
        ]),
        # dbc.Row([
        #     dbc.Col(html.P("The chosen equation is :", className="mb-0"),
        #             width="auto"),
        #     dbc.Col(html.Div(id='equation_select0'), width="auto"),
        #
        # ]),

        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P("dt (sec)", style={"height": "auto",
                                              "margin-bottom": "auto"}),
                    dcc.Input(id='dt-param-vivo2', type='number', value=3),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(["Suntag length (aa)",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Length of the SunTag in amino "
                                "acid.",
                                target="faq_prot_length"),
                    dcc.Input(id='suntag-length-param-vivo2', type='number',
                              value=800),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(["Protein length (aa)",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                            style={"height": "auto",
                                "margin-bottom": "auto"}),
                    dbc.Tooltip("Length of the protein in amino "
                                "acid.",
                                target="faq_prot_length"),
                    dcc.Input(id='prot-length-param-vivo2', type='number',
                              value=800),
                ]),
            ]),
        ]),
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P(["Missing point",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_missing_point",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("How many continuous missing point can "
                                "we recover. If it is too much the track will "
                                "not be analysed. For now, missing point are "
                                "created by calculated the mean value.",
                                target="faq_missing_point"),
                    dcc.Input(id='missing_point_param_vivo2',
                              type='number',
                              value=5),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(["Number of suntag",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Number of suntag repetition.",
                                target="faq_prot_length"),
                    dcc.Input(id='repetition-suntag-param-vivo2',
                              type='number',
                              value=32),
                ]),
            ]),

            dbc.Col([

            ]),

        ]),

        dbc.Row(
            [
                # dbc.Label(" "),
                dbc.Checklist(options=[{"label": "Force the analysis even if "
                                                 "time is not continuous.",
                                        "value": 0}],
                              id="switches_force_analysis2",
                              switch=False, ),
            ]),
        # Tick box to choose if we use the first dot for the analysis
        dbc.Row(
            [
                # dbc.Label(" "),
                dbc.Checklist(options=[{"label": "Include the first dot of "
                                                 "the autocorrelation "
                                                 "in the fit.",
                                        "value": 0}],
                              id="switches_first_dot2",
                              switch=False, ),
            ]),

        html.Br(),

        # Generate Button and Spinner Side by Side
        dbc.Row([
            # Analysis ONE track and ID of the track
            dbc.Col([
                html.H5("Display of one track",
                        style={"text-align": "center"}),
                dbc.Row([
                    html.Div([
                        html.P("id of the track to analyse",
                               style={"height": "auto",
                                      "margin-bottom": "auto"}),
                        dcc.Input(id='id_track2',
                                  type='number',
                                  value=0),
                    ]),
                ]),
                dbc.Row([
                    dbc.Col([
                        dbc.Button('Analyse and display ONE track',
                                   id='analyse_show_button2',
                                   className="mr-2",
                                   style={"width": "300px"}, ),
                    ], width="auto"),

                    dbc.Col([
                        dbc.Spinner(
                            children=[html.Div(id="loading_track_plot")],
                            size="sm", color="primary", type="border",
                            spinner_style={"margin-left": "10px"}
                        )
                    ], width="auto"),
                ], align="center", style={"margin-top": "10px"}),
                html.Div(id="loading_output_track"),
                html.Div(id='loading_output1'),
                html.Div(id='loading_output2'),
                html.Div(id='loading_output3'),
            ], style={'border': 'ridge',
                      'background-color': '#defffb'}),

            # Analysis ALL tracks and filename
            dbc.Col([
                html.H5("Analyse all tracks",
                        style={"text-align": "center"}),
                dbc.Row([
                    html.Div([
                        html.P("File name to save", style={"height": "auto",
                                                           "margin-bottom": "auto"}),
                        dcc.Input(id='save-results-name-vivo', type='text',
                                  value='datas_results'),
                    ]),
                ]),
                dbc.Row([
                    dbc.Col([
                        dbc.Button('Analyse ALL Tracks',
                                   id='start-analyze-btn-vivo',
                                   className="mr-2",
                                   style={"width": "300px"}, ),
                    ], width="auto"),

                    dbc.Col([
                        dbc.Spinner(
                            children=[html.Div(id="loading_analysis_vivo")],
                            size="sm", color="primary", type="border",
                            spinner_style={"margin-left": "10px"}
                        )
                    ], width="auto"),
                ], align="center", style={"margin-top": "10px",}),
                html.Div(id='analyze-output-vivo'),
                dcc.Download(id="download-csv"),
            ], style={'border': 'ridge',
                      'background-color': '#e6ffde'}),
        ]),

        html.Br(),
        html.Br(),

        dbc.Row([
            html.H5("Display of one track",
                    style={"text-align": "center"}),
            dcc.Graph(id='plot_results'),
        ]),

        html.Br(),
        html.Br(),

    )

