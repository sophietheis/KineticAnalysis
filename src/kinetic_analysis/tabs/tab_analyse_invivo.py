import numpy as np
import pandas as pd


from dash import html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .app_function import (upload_csv,)
from ..analysis.analysis_track import (single_track_analysis,
                                       validate_equation,
                                       fit_function,
                                       check_track_validity)


def layout():
    items = [dbc.DropdownMenuItem(
        "((t - x) / (c * t ** 2)) * heaviside((t - x), 0)"),
             dbc.DropdownMenuItem("ax+b")]
    return (
        # Analyse track in vivo tab
        # File selection and display
        dbc.Row([
            html.H4("Upload data",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
            dbc.Col([
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file to analyse "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_analyze_vivo',
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
                                    html.Div(id="loading_data_vivo")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-vivo'),
                ]),
            ], width=3),
            dbc.Col([
                html.Label("DataFrame Visualisation"),
                html.Div(id="table-container", children=[]),
            ], width=5)
        ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Choose columns name
        dbc.Row([
            html.H4("Confirm column name for the analysis",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
            dbc.Col([
                html.Div([
                    html.P("Track ID column",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_track',
                              type='text',
                              value="TRACK_ID",
                              style={'width': '200px'}),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(["Time column",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_time_col",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_time',
                              type='text',
                              value="FRAME",
                              style={'width': '200px'}),
                    dbc.Tooltip(
                        "It corresponds to a FRAME column. This will be "
                        "multiply by the dt time step.",
                        target="faq_time_col"),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("Intensity column",
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dcc.Input(id='col_intensity',
                              type='text',
                              value="MEAN_INTENSITY_CH1",
                              style={'width': '200px'}),
                ]),
            ]),
        ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Choose equation for the analysis
        dbc.Row([html.H4("Choose equation for the analysis",
                         style={"text-align": "center",
                                "color": "#10D79B"}),
                 html.Br(),
                 dcc.Markdown('''
                                 By default the model used for analyse track is :
                                 $$ \\frac{T - x}{cT^2}
                                 H(T - x), $$
                                 where $$c$$ is the initiation rate, and $$T$$ is the residence time.

                                 $$T=M/k$$ where $M$ is the RNA size (aa) and $k$ is the
                                 elongation rate.
                                ''',
                              mathjax=True),
                 html.Div([
                     html.P(["Equation ",
                             html.Span(className="fas fa-question-circle",
                                       id="faq_equation",
                                       style={"cursor": "pointer",
                                              "marginLeft": "5px"})],
                            style={"height": "auto",
                                   "margin-bottom": "auto"}),
                     dbc.Tooltip(
                         "You can change the equation. Symbols used are"
                         " x (fluorescence input), t (elongation rate) and "
                         "c (initiation rate). If you need other symbols, "
                         "please contact the admin.",
                         target="faq_equation"),

                     dbc.Row([
                         dbc.Col([
                             dbc.Select(["None",
                                         "((t - x) / (c * t ** 2)) * "
                                         "heaviside((t - x), 0)",
                                         "ax+b (not working"],
                                        "None",
                                        id="equation",
                                        className="mr-2",
                                        style={"width": 350,
                                               }, )
                         ], width=4),
                         dbc.Col([
                             dbc.Button('Valid equation',
                                        id='submit-button-equation',
                                        className="mr-2",
                                        style={"width": "150px"}, ),
                         ], width=3),
                     ]),
                     html.Br(),
                     dbc.Row([
                         dbc.Col([
                             dbc.Input(id='equation0',
                                       type='string',
                                       value="",
                                       placeholder="Enter equation",
                                       style={'width': '350px'}),
                             dbc.FormFeedback('Valid equation input!',
                                              type='valid'),
                             dbc.FormFeedback('Invalid equation input!',
                                              type='invalid'),
                         ], width=4),
                         dbc.Col([
                             dbc.Button('Valid equation',
                                        id='submit-button-equation0',
                                        className="mr-2",
                                        style={"width": "150px"}, ),
                         ], width=3),
                     ]),
                     # display the chosen equation and that it is valid
                     dbc.Row([
                         html.Div(id='loading_equation_output'),
                     ]),
                 ]),
                 ]),

        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Choose parameters for the analysis
        dbc.Row([
            html.H4("Confirm parameters for the analysis",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
        ]),

        # Put here the chosen equation
        dbc.Row([
            dbc.Col(html.P("The chosen equation is :", className="mb-0"),
                    width="auto"),
            dbc.Col(html.Div(id='equation_select'), width="auto"),
        ]),

        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P("dt (sec)", style={"height": "auto",
                                              "margin-bottom": "auto"}),
                    dcc.Input(id='dt-param-vivo', type='number', value=3),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(["Protein length (aa)",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Length of the protein + SunTag in amino "
                                "acid.",
                                target="faq_prot_length"),
                    dcc.Input(id='prot-length-param-vivo', type='number',
                              value=800),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P("File name to save", style={"height": "auto",
                                                       "margin-bottom": "auto"}),
                    dcc.Input(id='save-results-name-vivo', type='text',
                              value='datas_results'),
                ]),
            ]),
            html.Br(),
        ]),

        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P(["Missing point",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_missing_point",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("How many continuous missing point can "
                                "we recover. If it is too much the track will "
                                "not be analysed. For now, missing point are "
                                "created by calculated the mean value.",
                                target="faq_missing_point"),
                    dcc.Input(id='missing_point_param_vivo',
                              type='number',
                              value=5),
                ]),
            ]),
        ]),

        dbc.Row(
            [
                # dbc.Label(" "),
                dbc.Checklist(options=[{"label": "Force the analysis even if "
                                                 "time is not continuous.",
                                        "value": 0}],
                              id="switches_force_analysis",
                              switch=False, ),
            ]),
        # Tick box to choose if we use the first dot for the analysis
        dbc.Row(
            [
                # dbc.Label(" "),
                dbc.Checklist(options=[{"label": "Include the first dot of "
                                                 "the autocorrelation curve "
                                                 "in the fit.",
                                        "value": 0}],
                              id="switches_first_dot",
                              switch=False, ),
            ]),


        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        dbc.Row([
            html.Div([
                dcc.Checklist(
                    id='checkbox_simu',
                    options=[
                        {'label': 'Tracks from simulation', 'value':
                            'checked'}],
                    value=[],  # default is unchecked
                    labelStyle={'display': 'inline-block',
                                'margin-right': '10px'}
                )
            ])
        ]),
        dbc.Row([
            html.Br(),
            html.Br(),
        ]),
        # Generate Button and Spinner Side by Side
        dbc.Row([
            html.Br(),
            html.Br(),
            dbc.Col([
                dbc.Button('Start Analyse Tracks',
                           id='start-analyze-btn-vivo',
                           className="mr-2",
                           style={"width": "300px"}, ),
            ], width="auto"),

            dbc.Col([
                dbc.Spinner(
                    children=[html.Div(id="loading_analysis_vivo")],
                    size="sm",
                    color="primary",
                    type="border",
                    spinner_style={"margin-left": "10px"}
                )
            ], width="auto"),
            html.Div(id='analyze-output-vivo'),
            dcc.Download(id="download-csv")
        ])

    )


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-vivo', 'children'),
        Output("table-container", "children"),
        Output('loading_data_vivo', 'children'),
        Input('browse_directory_analyze_vivo', 'contents'),
    )
    def browse_directory_analyze_vivo(contents):
        df, output =  upload_csv(contents, app)
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
        Output("equation_select", "children", allow_duplicate=True),
        Input('submit-button-equation', 'n_clicks'),
        Input('equation', 'value'),
        prevent_initial_call=True,
    )
    def validate_input(n_clicks, value):
        if n_clicks:
            if value is None:
                return "None"
            else:
                v_bool, v_message, func_, expr_ = validate_equation(value)
                if value and v_bool:
                    app.data["equation_f"] = func_

                    return str(expr_)
                else:
                    return ""
        return ""

    @app.callback(
        Output('equation0', 'valid'),
        Output('equation0', 'invalid'),
        Output("loading_equation_output", "children"),
        Output("equation_select", "children"),
        Input('submit-button-equation0', 'n_clicks'),
        State('equation0', 'value')
    )
    def validate_input(n_clicks, value):
        if n_clicks:
            if not value:
                return False, False, f"Write an equation.", ""
            v_bool, v_message, func_, expr_ = validate_equation(value)
            if value and v_bool:
                app.data["equation_f"] = func_
                print(v_message)
                return True, False, v_message, str(expr_)
            else:
                return False, True, v_message, ""
        return False, False, "", ""

    @app.callback(
        Output('analyze-output-vivo', 'children'),
        Output('loading_analysis_vivo', 'children'),
        Output('download-csv', 'data'),
        Input('start-analyze-btn-vivo', 'n_clicks'),

        State('col_track', 'value'), #0
        State('col_time', 'value'), #1
        State('col_intensity', 'value'), #2
        State('dt-param-vivo', 'value'), #3
        State('prot-length-param-vivo', 'value'), #4
        State("missing_point_param_vivo", 'value'), #5
        State("switches_first_dot", "value"), #6
        State("switches_force_analysis", "value"), #7
        State('save-results-name-vivo', 'value'), #8
        State('checkbox_simu', 'value') #9
    )
    def start_analyze_tracks(n_clicks, *params):

        if n_clicks:
            print("click")
            if app.data['csv_to_analyse'] is None:
                return "No CSV file uploaded.", None, None

            try:
                print("start")
                # Read csv file
                df = app.data['csv_to_analyse']
                df.rename(columns={params[0]: 'TRACK_ID',
                                   params[1]: 'FRAME',
                                   params[2]: 'MEAN_INTENSITY_CH1',
                                  },
                         inplace=True)
                dt = float(params[3])
                prot_length = int(params[4])
                nb_missing_point = int(params[5])

                first_dot = bool(params[6])
                force_analysis = bool(params[7])
                # check_simu = 'checked' in params[9]

                ids_track = np.unique(df["TRACK_ID"])
                first_time = True
                # Analyse all tracks and save it
                for i in ids_track:
                    datas2 = df[(df["TRACK_ID"] == i)]

                    (valid,
                     x,
                     y,
                     x_fix,
                     y_fix) = check_track_validity(datas2,
                                                   i,
                                                   normalise_intensity=1,
                                                   delta_t=dt,
                                                   rtol=1e-1,
                                                   nb_missing_point=nb_missing_point,
                                                   )

                    if valid or force_analysis:
                        comment = ""
                        if force_analysis:
                            comment = "analysis forced"
                        (x_auto,
                         y_auto,
                         elongation_r,
                         translation_init_r,
                         perr) = single_track_analysis(x_fix,
                                                       y_fix,
                                                       delta_t=dt,
                                                       protein_size=prot_length,
                                                       mm=None,
                                                       normalise_auto=True,
                                                       method="original",
                                                       first_dot=first_dot,
                                                       simulation=False,
                                                       func_=app.data["equation_f"]
                                                       )
                        if first_time:
                            results = pd.DataFrame({"elongation_r": elongation_r,
                                                    "init_translation_r": translation_init_r,
                                                    "dt": dt,
                                                    "id": i,
                                                    "length": len(x_fix),
                                                    "perr0":perr[0],
                                                    "perr1": perr[1],
                                                    "comment": comment},
                                                   index=[0])
                            first_time = False

                        else:
                            results = pd.concat([results,
                                                 pd.DataFrame(
                                                     {"elongation_r": elongation_r,
                                                      "init_translation_r": translation_init_r,
                                                      "dt": dt,
                                                      "id": i,
                                                      "length": len(x_fix),
                                                      "perr0":perr[0],
                                                      "perr1": perr[1],
                                                      "comment":comment}, index=[0])
                                                 ], ignore_index=True)
                    else:
                        if first_time:
                            results = pd.DataFrame(
                                {"elongation_r": np.nan,
                                 "init_translation_r": np.nan,
                                 "dt": np.nan,
                                 "id": i,
                                 "length": len(x_fix),
                                 "perr0": np.nan,
                                 "perr1": np.nan,
                                 "comment":"cant be analysed"},
                                index=[0])
                            first_time = False

                        else:
                            results = pd.concat([results,
                                                 pd.DataFrame(
                                                     {
                                                         "elongation_r":
                                                             np.nan,
                                                         "init_translation_r": np.nan,
                                                         "dt": np.nan,
                                                         "id": i,
                                                         "length": len(x_fix),
                                                         "perr0": np.nan,
                                                         "perr1": np.nan,
                                                     "comment":"cant be "
                                                               "analysed"},
                                                     index=[0])
                                                 ], ignore_index=True)



                output_path = params[8] + ".csv"
                results.to_csv(output_path, index=False)

                return "Analysis completed and saved successfully!", None, dcc.send_file(output_path)
            except Exception as e:
                return f"Error: {str(e)}", None, None
        raise PreventUpdate


