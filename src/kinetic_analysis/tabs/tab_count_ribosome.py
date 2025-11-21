import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog

from dash import  html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_spinner

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from .app_function import (browse_directory)
from .app_function import (upload_csv)

from ..analysis.analyse_density import calculate_ribosome_density


def layout():
    return (html.Div([
        # Explanation at the begining of the page
        dbc.Row([
            html.P([
                "In this tab, you will estimate the number of ribosome on a "
                "translation site based on single protein fluorescence.",
                html.Br(),
                ]),
            dcc.Markdown('''
                                The calculation is based on this equation :
                                 $$ d = \\frac{F_{poly}}{F_{single}(L_{
                                 poi}+0.5L_{tag})}$$
                                 where $$F_{poly}$$ and $$F_{single}$$ are the 
                                 fluorescence of the polysome and single 
                                 protein respectively ; $$L_{poi}$$ and 
                                 $$L_{tag}$$ are the length of the protein 
                                 of interest and the tag respectively. 
                                ''',
                         mathjax=True),
        ]),
        html.Br(),
        # Upload dataframes
        html.H4("Upload data",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),
        html.H5("Single protein fluorescence",),
        dbc.Row([
            dbc.Col([
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file that contains single "
                               "protein fluorescence. "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_single_prot',
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
                                    html.Div(id="loading_data_single_prot")],
                                size="lm",  # "sm"
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-single_prot'),
                ]),
            ], width=3),
            dbc.Col([
                html.Label("Single protein fluorescence dataFrame visualisation"),
                html.Div(id="table_container_single", children=[]),
            ], width=5)
        ]),

        html.Br(),
        html.H5("Polysome fluorescence", ),
        dbc.Row([
            dbc.Col([
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file that contains polysome "
                               "fluorescence. "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_polysome',
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
                                    html.Div(id="loading_data_polysome")],
                                size="lm",  # "sm"
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-polysome'),
                ]),
            ], width=3),
            dbc.Col([
                html.Label(
                    "Polysome fluorescence dataFrame visualisation"),
                html.Div(id="table_container_polysome", children=[]),
            ], width=5)
        ]),

        html.Br(),
        html.Br(),
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
                    dcc.Input(id='param_prot_length_rib', type='number',
                              value=490,
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
                    dcc.Input(id='param_suntag_length_rib', type='number',
                              value=796,
                              style={'width': '200px'}),
                ]),


                html.Br(),
            ], width=3),

            # Calculate ribosome density
            dbc.Col([
                dbc.Button('Calculate ribosome density',
                           id='btn_calculate_ribosome_density',
                           className="mr-1"),
            ], width=5)
        ]),
        html.Div(id='output_ribosome_density'),
    ]),
    )


## Callbacks
def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-single_prot', 'children'),
        Output("table_container_single", "children"),
        Output('loading_data_single_prot', 'children'),
        Input('browse_directory_single_prot', 'contents'),
    )
    def browse_directory_single_prot(contents):
        df, output = upload_csv(contents, app, "csv_fluo_single")
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
        Output('selected-file-output-polysome', 'children'),
        Output("table_container_polysome", "children"),
        Output('loading_data_polysome', 'children'),
        Input('browse_directory_polysome', 'contents'),
    )
    def browse_directory_polysome(contents):
        df, output = upload_csv(contents, app, "csv_fluo_polysome")
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
        Output('output_ribosome_density', 'children'),
        Input('btn_calculate_ribosome_density', 'n_clicks'),
        State('param_prot_length_rib', 'value'),  #0
        State('param_suntag_length_rib', 'value'),  #1
    )
    def calculate_density(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """
        if n_clicks:
            print("click")
            try:
                L_poi = float(params[0])
                L_tag = float(params[1])

                m_intensity_single, result = calculate_ribosome_density(
                    app.data["csv_fluo_single"],
                    app.data["csv_fluo_polysome"],
                    L_poi,
                    L_tag)
                print("calculaiton done")
                out_string = ("Mean single protein intensity : " + str(
                    np.round(m_intensity_single,2)) +
                              "\n Mean polysome intensity : " + str(
                            np.round(result["INTENSITY"].mean(), 2)) +
                              "\n Mean ribosome density : " + str(
                    np.round(result["ribosome_density"].mean(),2)))
                print(out_string)
                return out_string
            except Exception as e:
                print(e)
                return "Problem"
        raise PreventUpdate
