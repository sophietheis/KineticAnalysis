from dash import html, dcc
import dash_bootstrap_components as dbc


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([
                "In this tab, you will estimate the "
                " optimal microscopy acquisition settings for a live-cell SunTag "
                "imaging experiment before data collection. ",
                html.Br(),
                "It is recommended that track length is at least 3 times "
                "the time to translate one (Suntag+protein).",
                html.Br(),
                "The recommended frame rate should be lower that 1/3 "
                "of the time to translate the SunTag ",
            ]),
        ]),

        html.Br(),

        dbc.Row([
            dbc.Col([
                # Select value
                html.H4(children="Input value",
                        style={"text-align": "center",
                               "color": "#10D79B"}),
            ], width=5),
            dbc.Col([
                # Select value
                html.H4(children="Output",
                        style={"text-align": "center",
                               "color": "#10D79B"}),
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
                dbc.Button(children='Calculate acquisition parameters',
                           id='btn_calculate_acquisition',
                           className="mr-1"),
            ], width=5),

            dbc.Col([
                html.Div(id='output_acquisition'),
            ])
        ]),

    ]),
    )
