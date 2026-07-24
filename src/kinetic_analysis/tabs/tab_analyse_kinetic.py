from dash import html, dcc
import dash_bootstrap_components as dbc


def layout():
    items_solver = ["",
                    "Exact equation",
                    "Approximate equation",
                    "Approximate epitope",
                    # "linear fit",
                    ]

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
                                id='browse_directory_analyze_vivo2',
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
                                    html.Div(id="loading_data_vivo2")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-vivo2'),
                ]),
            ], width=3),

            dbc.Col(children=[
                html.Label("DataFrame Visualisation"),
                html.Div(id="table-container2", children=[]),
            ], width=9)
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
                    html.P(children="Track ID column",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_track2',
                              type='text',
                              value="TRACK_ID",
                              style={'width': '200px'}),
                ]),
            ]),

            # Time column name
            dbc.Col([
                html.Div([
                    html.P(children=["Time column",
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
                        children="It corresponds to a FRAME column. "
                                 "This will be multiply by the dt time step.",
                        target="faq_time_col"),
                ]),
            ]),

            # Intensity column name
            dbc.Col([
                html.Div([
                    html.P(children="Intensity column",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_intensity2',
                              type='text',
                              value="MEAN_INTENSITY_CH1",
                              style={'width': '200px'}),
                ]),
            ]),
        ]),

        html.Br(),
        html.Br(),

        # CHOOSE THE EQUATION FOR THE ANALYSIS
        html.H4(children="Choose equation for the analysis",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),

        dbc.Row([
            dbc.Col([
                dbc.Select(options=items_solver,
                           id="choose-solver1",
                           className="mr-2",
                           style={"width": 250}, ),

                html.Br(),
                # Put here the chosen equation
                html.P(children="The chosen solver is: ",
                       className="mb-0"),
                html.Div(id="choosen-solver1"),

            ], width=3),

            dbc.Col([
                dcc.Markdown(children='''
                            If you choose **exact equation**, the default the 
                             equation used is :
                             $$
                             G(T)_{exact} = \\frac{Pe^{-kT}}{P^2(
                             \\frac{N(N+1)}{2}+NM)^2 }[ \sum_{
                             n=0}^{N}(N-n)(N-n+1)\\frac{(2N+n+1)}{6}\\frac{(
                             kT)^n}{n!}
                             $$
                             $$ 
                             + N\sum_{n=0}^{N}n\\frac{(2N-n+1}{
                             2}\\frac{(kT)^n}{n!}+ N^2\\frac{N+1}{2}\sum_{
                             n=N}^{M}\\frac{(kT)^n}{n!} 
                             $$
                             $$+ N\sum_{n=1}^{
                             N}n\\frac{(1+n)}{2}\\frac{(kT)^{M+N-n}}{(
                             M+N-n)!}+ N^2\sum_{n=0}^{M}(M-n)\\frac{(
                             kT)^n}{n!} ] 
                             $$
                             where $$c$$ is the initiation rate, $$k$$ 
                             the elongation rate and $$T$$ is 
                             the residence time. The steady state 
                             occupancy $$P$$ is defined by 
                             $$\\frac{c}{k}$$.
                             ''',
                             mathjax=True),
                html.Hr(style={'borderWidth': "0.3vh", "width": "25%",
                               "color": "#10D79B"}),
                dcc.Markdown(children='''
                             If you choose **approximate equation**, 
                             the default 
                             the equation used is :
                             $$ G(x) = \\frac{T - x}{cT^2}
                             H(T - x), $$
                             where $$c$$ is the initiation rate, and $$T$$ is 
                             the residence time.

                             $$T=M/k$$ where $M$ is the RNA size (aa) and $k$ 
                             is the elongation rate.
                             ''',
                             mathjax=True),
                html.Hr(style={'borderWidth': "0.3vh", "width": "25%",
                               "color": "#10D79B"}),
                dcc.Markdown(children='''             
                             If you choose **approximate epitope**, 
                             the default equation used is : 
                             $$ G(T) = \\frac{k}{c}(\\frac{2}{3})\\frac{
                             1}{(N(N+1))^2}e^{-kT} \\sum_{n=0}^N(N-n)(
                             N-n+1)(2N+n+1)\\frac{(kT)^n}{n!} $$
                            ''',
                             mathjax=True),
            ], width=9),
        ]),



        html.Br(),
        html.Br(),

        # CONFIRM PARAMETER FOR THE ANALYSIS
        html.H4(children="Confirm parameters for the analysis",
                style={"text-align": "center",
                       "color": "#10D79B"}),
        html.Br(),

        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P(children="dt (sec)",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='dt-param-vivo2', type='number', value=3),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(children=["Suntag length (aa)",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_suntag_length",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Length of the SunTag in amino "
                                "acid.",
                                target="faq_suntag_length"),
                    dcc.Input(id='suntag-length-param-vivo2', type='number',
                              value=800),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(children=["Protein length (aa)",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_prot_length",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Length of the protein in amino "
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
                    html.P(children=["Missing point",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_missing_point",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="How many continuous missing point can "
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
                    html.P(children=["Number of suntag",
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_suntag_number",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Number of suntag repetition.",
                                target="faq_suntag_number"),
                    dcc.Input(id='repetition-suntag-param-vivo2',
                              type='number',
                              value=32),
                ]),
            ]),

            dbc.Col([
            ]),

        ]),
        html.Br(),
        dbc.Row([
            dbc.Checklist(options=[{"label": "Force the analysis even if "
                                                "time is not continuous.",
                                    "value": 0}],
                            id="switches_force_analysis2",
                            switch=False, ),
        ]),
        html.Br(),            
        dbc.Row([
            dbc.Col([
                dbc.Checklist(options=[{"label": "Correct queuing",
                                                        "value": 0}],
                                              id="switches_correct_queuing_analysis2",
                                              switch=False, ),
            ]),
            dbc.Col([
                html.Div([
                    html.P(children=["Mean ribosome occupancy",
                                        html.Span(className="fas fa-question-circle",
                                                id="faq_rib_occupancy",
                                                style={"cursor": "pointer",
                                                        "marginLeft": "5px"})],
                            style={"height": "auto",
                                    "margin-bottom": "auto"}),
                    dbc.Tooltip(children="mean ribosome occupancy.",
                                target="faq_rib_occupancy"),
                    dcc.Input(id='rib_occupancy-param-vivo2',
                                type='number',
                                value=10),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(children=["Ribosome footprint",
                                        html.Span(className="fas fa-question-circle",
                                                id="faq_rib_footprint",
                                                style={"cursor": "pointer",
                                                        "marginLeft": "5px"})],
                            style={"height": "auto",
                                    "margin-bottom": "auto"}),
                    dbc.Tooltip(children="Ribosome footprint.",
                                target="faq_rib_footprint"),
                    dcc.Input(id='rib_footprint-param-vivo2',
                                type='number',
                                value=24),
                ]),
            ])
        ]),
        
        html.Br(),

        # Generate Button and Spinner Side by Side
        dbc.Row([
            # Analysis ONE track and ID of the track
            dbc.Col(children=[
                html.H5(children="Display of one track",
                        style={"text-align": "center"}),
                dbc.Row([
                    html.Div([
                        html.P(children="id of the track to analyse",
                               style={"height": "auto",
                                      "margin-bottom": "auto"}),
                        dcc.Input(id='id_track2',
                                  type='number',
                                  value=0),
                    ]),
                ]),
                dbc.Row(children=[
                    dbc.Col(children=[
                        dbc.Button(children='Analyse and display ONE track',
                                   id='analyse_show_button2',
                                   className="mr-2",
                                   style={"width": "300px"}, ),
                        dbc.Alert(
                            "Analysis did not work",
                            id="alert-analysis",
                            dismissable=True,
                            fade=False,
                            is_open=False,
                            color="danger",
                        ),
                    ], width="auto"),

                    dbc.Col(children=[
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
            dbc.Col(children=[
                html.H5(children="Analyse all tracks",
                        style={"text-align": "center"}),
                dbc.Row([
                    html.Div([
                        html.P(children="File name to save",
                               style={"height": "auto",
                                      "margin-bottom": "auto"}),
                        dcc.Input(id='save-results-name-vivo', type='text',
                                  value='datas_results'),
                    ]),
                ]),
                dbc.Row(children=[
                    dbc.Col(children=[
                        dbc.Button(children='Analyse ALL Tracks',
                                   id='start-analyze-btn-vivo',
                                   className="mr-2",
                                   style={"width": "300px"}, ),
                    ], width="auto"),

                    dbc.Col(children=[
                        dbc.Spinner(
                            children=[html.Div(id="loading_analysis_vivo")],
                            size="sm", color="primary", type="border",
                            spinner_style={"margin-left": "10px"}
                        )
                    ], width="auto"),
                ], align="center", style={"margin-top": "10px"}),
                html.Div(id='analyze-output-vivo'),
                dcc.Download(id="download-csv"),
            ], style={'border': 'ridge',
                      'background-color': '#e6ffde'}),
        ]),

        html.Br(),
        html.Br(),

        # DISPLAY RESULTS OF ONE TRACK
        dbc.Row([
            html.H5(children="Display of one track curves and result",
                    style={"text-align": "center"}),
            dcc.Graph(id='plot_results'),
        ]),


    )
