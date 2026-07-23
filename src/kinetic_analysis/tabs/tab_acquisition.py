from dash import html, dcc
import dash_bootstrap_components as dbc


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([
                "In this tab, you will calculate the ideal acquisition parameters for your experiment. ",
                html.Br(),
            ]),
        ]),
        
        html.Br(),
        # Select value
        html.H4(children="Confirm value",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),
        dbc.Row([
            # Track ID column name
            dbc.Col([
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
            ]),

            # X column name
            dbc.Col([
                html.Div([
                    html.P(children="Suntag length (aa) ",
                            style={"height": "auto",
                                    "margin-bottom": "auto"}),
                    dcc.Input(id='param_suntag_length_acquisition',
                                type='number',
                                value=796,
                                style={'width': '200px'}),
                ]),
            ]),

            # Y column name
            dbc.Col([
                html.Div([
                    html.P(children="Translation rate (aa/sec) ",
                            style={"height": "auto",
                                    "margin-bottom": "auto"}),
                    
                    dcc.Input(id='param_translation_rate_acquisition',
                                type='number',
                                value=24,
                                style={'width': '200px'}),
                ]),
            ]),

            # Z column name
            dbc.Col([
                 html.Div([
                    html.P(children="Number of full translated protein",
                            style={"height": "auto",
                                    "margin-bottom": "auto"}),
                   
                    dcc.Input(id='param_full_translated_protein_acquisition',
                                type='number',
                                value=1,
                                style={'width': '200px'}),
                ]),
            ]),

            # A column name
            dbc.Col([
                    html.Div([
                    html.P(children="Number of stem loops",
                            style={"height": "auto",
                                    "margin-bottom": "auto"}),
                    dcc.Input(id='param_stem_loops_acquisition',
                                type='number',
                                value=24,
                                style={'width': '200px'}),
                ]),
            ]),
        ]),

        html.Br(),
        html.Br(),
        dbc.Row([
            dbc.Col([
                dbc.Button(children='Calculate acquisition parameters',
                           id='btn_calculate_acquisition',
                           className="mr-1"),
            ], width=5)
        ]),
        html.Div(id='output_acquisition'),
    ]),
    )    
