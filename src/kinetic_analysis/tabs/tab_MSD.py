from dash import html, dcc
import dash_bootstrap_components as dbc


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([
                "In this tab, you will calculate the MSD (Mean Square "
                "Displacement) of fluorescence dots in 3D.",
                html.Br(),
                ]),
            dcc.Markdown('''
                                The calculation is based on this equation :
                                 $$MSD = <\sum{(r_{t}-r_{t-1})²}>$$
                                 where $$r_{t}$$ and $$r_{t-1}$$ are the 
                                 position at time t and t-1. 
                                ''',
                         mathjax=True),
            html.P([
                "Columns name should be \"x\", \"y\", \"z\"",
                html.Br(),
            ]),
        ]),
        html.Br(),
        # Upload dataframes
        html.H4("Upload data",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),

        dbc.Row([
            dbc.Col([
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file. "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_msd',
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
                                    html.Div(id="loading_data_msd")],
                                size="lm",  # "sm"
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-msd'),
                ]),

                html.Br(),
                html.Div(id="table_container_msd", children=[]),
                html.Br(),
            ], width=6),

        ]),
        html.Br(),
        html.Br(),
        # CHOOSE COLUMN NAME
        html.H4(children="Confirm column name for the analysis",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),
        dbc.Row([
            # Track ID column name
            dbc.Col([
                html.Div([
                    html.P(children="Track ID",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_ID',
                              type='text',
                              value="TRACK_ID",
                              style={'width': '200px'}),
                ]),
            ]),

            # X column name
            dbc.Col([
                html.Div([
                    html.P(children="x",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_x',
                              type='text',
                              value="POSITION_X",
                              style={'width': '200px'}),
                ]),
            ]),

            # Y column name
            dbc.Col([
                html.Div([
                    html.P(children="y",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_y',
                              type='text',
                              value="POSITION_Y",
                              style={'width': '200px'}),
                ]),
            ]),

            # Z column name
            dbc.Col([
                html.Div([
                    html.P(children="z",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_z',
                              type='text',
                              value="POSITION_Z",
                              style={'width': '200px'}),
                ]),
            ]),
        ]),

        html.Br(),
        html.Br(),
        dbc.Row([
            # Calculate MSD
            dbc.Col([
                dbc.Button('Calculate MSD',
                           id='btn_calculate_MSD',
                           className="mr-1"),
            ], width=5)
        ]),
        html.Div(id='output_MSD'),
        dcc.Download(id="download_csv_msd"),
    ]),
    )
