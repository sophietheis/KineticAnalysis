from dash import html, dcc
import dash_bootstrap_components as dbc

from kineticanalysis.utils.texts import t
from kineticanalysis.tabs.tab_utils import (button,
                                            input_case_with_tooltip,
                                            input_case,
                                            title_h4,
                                            color_line)


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
        title_h4("Upload data"),
        html.Br(),
        dbc.Row([
            dbc.Col(children=[
                html.Div([
                    # Choose a csv file
                    html.Label(t("tab_analysis.select_csv")),
                    html.Br(),
                    dbc.Row([
                        dbc.Col(children=[
                            dcc.Upload(
                                id='browse_directory_analyze_vivo2',
                                children=button(t("button.upload_csv")),
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
        title_h4("Confirm column name for the analysis"),
        html.Br(),
        dbc.Row([
            # Track ID column name
            dbc.Col([
                input_case(child=t("labels.Track_id"),
                           input_id='col_track_kinetic',
                           input_type='text',
                           input_value=t("default_input_value.id")),
            ]),

            # Time column name
            dbc.Col([
                input_case_with_tooltip(child=t("labels.Time"),
                                        input_id='col_time_kinetic',
                                        tooltip_id="faq_time_col_analysis",
                                        tooltip_text=t("tooltips.time_col"),
                                        input_type='text',
                                        input_value=t("default_input_value.time")),
            ]),

            # Intensity column name
            dbc.Col([
                input_case(child=t("labels.intensity"),
                           input_id='col_intensity_kinetic',
                           input_type='text',
                           input_value=t("default_input_value.intensity")),
            ]),
        ]),

        html.Br(),
        html.Br(),

        # CHOOSE THE EQUATION FOR THE ANALYSIS
        title_h4("Choose equation for the analysis"),

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
                html.P(children="If you choose **exact equation**, the default the equation used is :"),
                dcc.Markdown(children=t("equations.g_exact"), mathjax=True),
                color_line(),
                html.P(children="If you choose **approximate equation**, the default the equation used is :"),
                dcc.Markdown(children=t("equations.g_approx"), mathjax=True),
                color_line(),
                html.P(children="If you choose **approximate epitope**, the default the equation used is :"),
                dcc.Markdown(children=t("equations.g_epitope"), mathjax=True),
            ], width=9),
        ]),

        html.Br(),
        html.Br(),

        # CONFIRM PARAMETER FOR THE ANALYSIS
        title_h4("Confirm parameters for the analysis"),
        html.Br(),

        dbc.Row([
            dbc.Col([
                input_case(child=t("labels.dt"),
                           input_id='dt-param-vivo_kinetic',
                           input_type='number',
                           input_value=t("default_input_value.dt_analysis")),
            ]),
            dbc.Col([
                input_case_with_tooltip(child=t("labels.suntag"),
                                        input_id='suntag-length-param-vivo_kinetic',
                                        tooltip_id="faq_suntag_length",
                                        tooltip_text=t("tooltips.suntag_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.suntag_length")),
            ]),
            dbc.Col([
                input_case_with_tooltip(child=t("labels.prot"),
                                        input_id='prot-length-param-vivo_kinetic',
                                        tooltip_id="faq_prot_length",
                                        tooltip_text=t("tooltips.prot_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.prot_length")),
            ]),
        ]),
        dbc.Row([
            dbc.Col([
                html.Div([
                    input_case_with_tooltip(child=t("labels.missing_point"),
                                            input_id='missing_point_param_vivo_kinetic',
                                            tooltip_id="faq_missing_point",
                                            tooltip_text=t("tooltips.missing_point"),
                                            input_type='number',
                                            input_value=t("default_input_value.missing_point")),
                ]),
            ]),
            dbc.Col([
                input_case_with_tooltip(child=t("labels.nb_stemloop"),
                                        input_id='repetition-suntag-param-vivo_kinetic',
                                        tooltip_id="faq_suntag_number",
                                        tooltip_text=t("tooltips.nb_stemloop"),
                                        input_type='number',
                                        input_value=t("default_input_value.nb_stemloop")),
            ]),

            dbc.Col([
            ]),

        ]),
        html.Br(),
        dbc.Row([
            dbc.Checklist(options=[{"label": "Force the analysis even if "
                                             "time is not continuous.",
                                    "value": 0}],
                          id="switches_force_analysis_kinetic",
                          switch=False, ),
        ]),
        html.Br(),
        dbc.Row([
            dbc.Col([
                dbc.Checklist(options=[{"label": "Correct queuing",
                                                 "value": 0}],
                              id="switches_correct_queuing_analysis_kinetic",
                              switch=False, ),
            ]),
            dbc.Col([
                input_case_with_tooltip(child=t("labels.mean_rib_occ"),
                                        input_id='rib_occupancy-param-vivo_kinetic',
                                        tooltip_id="faq_rib_occupancy",
                                        tooltip_text=t("tooltips.mean_rib_occ"),
                                        input_type='number',
                                        input_value=t("default_input_value.mean_rib_occ")),
            ]),
            dbc.Col([
                input_case_with_tooltip(child=t("labels.footprint"),
                                        input_id='rib_footprint-param-vivo_kinetic',
                                        tooltip_id="faq_rib_footprint",
                                        tooltip_text=t("tooltips.footprint"),
                                        input_type='number',
                                        input_value=t("default_input_value.footprint_analysis")),
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
                    input_case(child=t("labels.single_track_id"),
                               input_id='id_track_kinetic',
                               input_type='number',
                               input_value=t("default_input_value.track_id")),
                ]),
                dbc.Row(children=[
                    dbc.Col(children=[
                        button(child='Analyse and display ONE track',
                               id_='analyse_show_button_kinetic',
                               width="300px"),
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
                    input_case(child=t("labels.filename"),
                               input_id='save-results-name-vivo',
                               input_type='number',
                               input_value=t("default_input_value.filename_results")),
                ]),
                dbc.Row(children=[
                    dbc.Col(children=[
                        button(child='Analyse ALL Tracks',
                               id_='start-analyze-btn-vivo',
                               width="300px"),
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
