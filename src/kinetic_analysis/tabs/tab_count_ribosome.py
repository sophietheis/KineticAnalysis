from dash import html, dcc
import dash_bootstrap_components as dbc


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
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

        dbc.Row([
            dbc.Col([
                html.H5("Single protein fluorescence", ),
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
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-single_prot'),
                ]),

                html.Br(),
                html.Label("Single protein fluorescence dataFrame visualisation"),
                html.Div(id="table_container_single", children=[]),
                html.Br(),
            ], width=6),

            dbc.Col([
                html.H5("Polysome fluorescence", ),
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

                html.Br(),
                html.Label(
                    "Polysome fluorescence dataFrame visualisation"),
                html.Div(id="table_container_polysome", children=[]),

            ], width=6),
        ]),
        html.Br(),
        html.Br(),
        dbc.Row([
            html.Label("After uploading your data, select a column in each "
                       "dataframe that corresponds to the intensity."),
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
        dcc.Download(id="download_csv"),
    ]),
    )
