from dash import html, dcc
import dash_bootstrap_components as dbc

from kineticanalysis.tabs.tab_utils import title_h4, input_case, button
from kineticanalysis.utils.texts import t


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([t("tab_acquisition.intro"),
                    html.Br(),
                    t("tab_acquisition.intro1"),
                    html.Br(),
                    t("tab_acquisition.intro2"),
                    ]),
        ]),
        html.Br(),

        dbc.Row([
            dbc.Col([
                title_h4("Input value"),
            ], width=5),
            dbc.Col([
                title_h4("Output"),
            ], width=5),
        ]),
        html.Br(),

        dbc.Row([
            dbc.Col([
                # Track ID column name
                html.Div([
                    html.P(children="Protein length (aa) ",
                           style={"height": "auto",
                                    "margin-bottom": "auto"}
                           ),
                    dcc.Input(id='param_prot_length_acquisition',
                              type='number',
                              value=490,
                              style={'width': '200px'}),
                ]),

                # X column name
                html.Div([
                    html.P(children="Suntag length (aa) ",
                           style={"height": "auto",
                                    "margin-bottom": "auto"}),
                    dcc.Input(id='param_suntag_length_acquisition',
                              type='number',
                              value=796,
                              style={'width': '200px'}),
                ]),

                # A column name
                html.Div([
                    html.P(children="Number of stem loops",
                           style={"height": "auto",
                                    "margin-bottom": "auto"}),
                    dcc.Input(id='param_stem_loops_acquisition',
                              type='number',
                              value=24,
                              style={'width': '200px'}),
                ]),

                html.Br(),

                # Y column name
                html.Div([
                    html.P(children="Translation rate (aa/sec) ",
                           style={"height": "auto",
                                    "margin-bottom": "auto"}),

                    dcc.Input(id='param_translation_rate_acquisition',
                              type='number',
                              value=24,
                              style={'width': '200px'}),
                ]),

                # Z column name
                html.Div([
                    html.P(children="Number of full translated protein",
                           style={"height": "auto",
                                    "margin-bottom": "auto"}),

                    dcc.Input(id='param_full_translated_protein_acquisition',
                              type='number',
                              value=1,
                              style={'width': '200px'}),
                ]),
                html.Br(),
                button(child=t("button.acquisition"),
                       id_='btn_calculate_acquisition'),
            ], width=5),

            dbc.Col([
                html.Div(id='output_acquisition'),
            ])
        ]),

    ]),
    )
