from dash import html, dcc
import dash_bootstrap_components as dbc


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.Br(),
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

        # Define parameter of the simulation
        dbc.Row([
            dbc.Col(children=[
                html.Div([
                    html.P(children=["Protein length (aa) ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_prot_length",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}
                           ),
                    dbc.Tooltip(children="Length of the protein in amino "
                                         "acid. This value will be added "
                                         "to suntag length.",
                                target="faq_param_prot_length"),

                    dcc.Input(id='param_prot_length',
                              type='number',
                              value=490,
                              style={'width': '200px'}),
                ]),

                html.Div([
                    html.P(children=["Suntag length (aa) ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_suntag_length",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Length of the suntag in amino "
                                         "acid. This value will be added "
                                         "to protein length.",
                                target="faq_param_suntag_length"),
                    dcc.Input(id='param_suntag_length',
                              type='number',
                              value=796,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Number of suntag ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_nb_suntag",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Number of suntag repetition",
                                target="faq_param_nb_suntag"),
                    dcc.Input(id='param_nb_suntag',
                              type='number',
                              value=32,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Fluorescence one suntag ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_fluo_one_suntag",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Fluorescence of one suntag, use as "
                                         "reference for fluorescence profile.",
                                target="faq_param_fluo_one_suntag"),
                    dcc.Input(id='param_fluo_one_suntag',
                              type='number',
                              value=4,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Translation rate (aa/sec) ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_translation_rate",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Translation rate.",
                                target="faq_param_translation_rate"),
                    dcc.Input(id='param_translation_rate',
                              type='number',
                              value=24,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Initiation rate (ribosome/sec) ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_initiation_rate",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={
                               "height": "auto",
                               "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Initiation rate.",
                                target="faq_param_initiation_rate"),
                    dcc.Input(id='param_initiation_rate',
                              type='number',
                              value=1,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Retention time (sec) ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_retention_time",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Protein stay at the RNA for some time",
                                target="faq_param_retention_time"),
                    dcc.Input(id='param_retention_time',
                              type='number',
                              value=0,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Suntag position (begin or end) ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_pos_suntag",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Choose if Suntag is before or "
                                         "after the protein.",
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
                    html.P(children=["Noise",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_noise",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Add random noise to the protein "
                                         "fluorescence profile. "
                                         "0 if no noise. ",
                                target="faq_param_noise"),
                    dcc.Input(id='param_noise', type='number', value=0.,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["dt (sec)",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_dt",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(
                        children="Time step use to generate track. "
                                 "This value should be quite small "
                                 "regarding the time step used for "
                                 "the analysis",
                        target="faq_param_dt"),
                    dcc.Input(id='param_dt',
                              type='number',
                              value=0.1,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Length of one track (sec)",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_length",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(
                        children="Length duration of one track."
                        "used for the analysis",
                        target="faq_param_length"),
                    dcc.Input(id='param_length',
                              type='number',
                              value=6000,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children=["Number of tracks ",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_nb_tracks",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Number of tracks to be generated.",
                                target="faq_param_nb_tracks"),
                    dcc.Input(id='param_nb_tracks',
                              type='number',
                              value=100,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(children="File name to save",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='param_filename',
                              type='text',
                              value='datas',
                              style={'width': '200px'}),
                ]),
            ], width=4),

            # Graph area
            dbc.Col(children=[
                dcc.Graph(id='profile-plot'),
            ], width=8),
        ]),

        html.Br(),
        html.Br(),

        # Button row
        dbc.Row([
            # Show plot profile button
            dbc.Col(children=[
                dbc.Button(children='Show Profile',
                           id='show-profile-btn',
                           className="mr-1"),
            ], width=2),

            # Generate track Button and Spinner Side by Side
            dbc.Col([
                html.Div([
                    dbc.Col(children=[
                        dcc.Store(id="start", data=""),
                        dcc.Store(id="complete", data=""),
                        dbc.Button(dbc.Spinner(
                            html.Span(children="Generate tracks",
                                      id="loading_generate")),
                            id="start-gen-tracks-btn"),
                    ], width=3),
                ]),
                html.Div(id='gen-tracks-output'),
                dcc.Download(id="download-csv2"),
            ]),
        ]),
        # dbc.Row([dcc.Interval(id="progress-generate",
        #                       n_intervals=0,
        #                       interval=1000,)]),
    ]),
    )
