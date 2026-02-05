from dash import html, dcc
import dash_bootstrap_components as dbc


def layout():

    return (
        html.Br(),

        # FILE IMPORT AND DISPLAY
        html.H4(children="Upload data",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),
        dbc.Row([
            dbc.Col(children=[
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file to analyse "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col(children=[
                            dcc.Upload(
                                id='browse_directory_combine',
                                children=dbc.Button(children='Upload csv file',
                                                    className="mr-2",
                                                    style={"width": "150px"}
                                                    ),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col(children=[
                            dbc.Spinner(
                                children=[
                                    html.Div(id="loading_data_combine")],
                                size="lm",  #"sm"
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-combine'),
                ]),
            ], width=3),

            dbc.Col(children=[
                html.Label("DataFrame Visualisation"),
                html.Div(id="table-container_combine", children=[]),
            ], width=9)
        ]),

        html.Br(),

        dbc.Row([
            html.P(children=["There is : ",
                             html.Span(id='nb_tracks_init', children=''),
                             " track(s)."
        ],
                   className="mb-0"),
        ]),

        html.Br(),

        # # input values
        # html.H4(children="Confirm column name for the analysis",
        #         style={"text-align": "center",
        #                "color": "#10D79B"}),
        # html.Br(),
        dbc.Row([
            # nb track in the combined track
            dbc.Col([
                html.Div([
                    html.P(children="How many tracks in the combined tracks?",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='nb_tracks',
                              type='number',
                              value="2",
                              style={'width': '200px'}),
                ]),
            ]),

            # number of new tracks created
            dbc.Col([
                html.Div([
                    html.P(children="How many new track created ?",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='nb_new_tracks',
                              type='number',
                              value="10",
                              style={'width': '200px'}),
                ]),
            ]),
        ]),

        html.Br(),
        dbc.Row([
            dbc.Col([
                html.Div([
                    dbc.Col(children=[
                        dcc.Store(id="start2", data=""),
                        dcc.Store(id="complete2", data=""),
                        dbc.Button(dbc.Spinner(
                            html.Span(children="Combine tracks",
                                      id="loading_combination")),
                            id="combine-tracks-btn"),
                    ], width=3),
                ]),
                html.Div(id='combine-tracks-output'),
                dcc.Download(id="download-csv3"),
            ]),
        ]),

        html.Br(),

    )
