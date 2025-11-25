import numpy as np

from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc


import plotly.graph_objs as go
from plotly.subplots import make_subplots


from ..analysis.contribution import calculate_contribution


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([
                "In this tab, you will be able to visualise the "
                "contribution of each term in the equation in order "
                "to help you choose the best method of analysis. ",
                html.Br(),
                "This is based on the work of Larson et al. ",
                html.Br(),
                ]),
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
                    html.P(["Estimate elongation rate (aa/sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_translation_rate",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Elongation rate.",
                                target="faq_param_translation_rate"),
                    dcc.Input(id='param_elongation_rate', type='number',
                              value=0.3,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Estimate initiation rate (ribosome/sec) ",
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
                    dcc.Input(id='param_initiation_rate', type='number',
                              value=0.015,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["tau",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_dt",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(
                        "Maximum length",
                        target="faq_param_dt"),
                    dcc.Input(id='param_tau', type='number', value=150,
                              style={'width': '200px'}),
                ]),
                html.Br(),
            ], width=4),

            # Show plot contribution
            dbc.Col([
                dcc.Graph(id='equation-plot'),

            ], width=8)
        ]),
        dbc.Row([
            dbc.Button('Show Contribution', id='show-contribution-btn',
                       className="mr-1"),
        ]),
    ]),
    )

