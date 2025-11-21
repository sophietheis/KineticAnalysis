import numpy as np
import pandas as pd

from dash import  html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from .app_function import (upload_csv)
from ..analysis.analysis_track import (single_track_analysis,
                                       validate_equation,
                                       fit_function,
                                       fit_function_exact,
                                       check_track_validity)


def layout():
    items_equation = ["",
                      "((t - x) / (c * t ** 2)) * heaviside((t - x), 0)",
                      ]
    items_solver = ["",
                    "exact_fit",
                    "approximation fit",
                    "linear fit"]

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
                                id='browse_directory_analyze_vivo2',
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
                                    html.Div(id="loading_data_vivo2")],
                                size="lm", #"sm"
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-vivo2'),
                ]),
            ], width=3),
            dbc.Col([
                html.Label("DataFrame Visualisation"),
                html.Div(id="table-container2", children=[]),
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
                    dcc.Input(id='col_track2',
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
                    dcc.Input(id='col_time2',
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
                    dcc.Input(id='col_intensity2',
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
        dbc.Row([
            html.H4("Choose equation and solver for the analysis",
                    style={"text-align": "center",
                           "color": "#10D79B"}),
            html.Br(),
            # Solver
            dbc.Col([
                html.H5("Solver"),
                dcc.Markdown('''
                     Choose the solver. \n
                     You can choose between exact_fit, approximation_fit 
                     and linear_fit. \n 
                    ''',
                             mathjax=True),

                dbc.Select(options=items_solver,
                           id="choose-solver1",
                           className="mr-2",
                           style={"width":350},),

                dcc.Markdown('''
                            If you choose exact_fit, the default the 
                             equation used is :
                             $$ G(T) = \\frac{k}{c}(\\frac{2}{3})\\frac{
                             1}{(N(N+1))^2}e^{-kT} \\sum_{n=0}^N(N-n)(
                             N-n+1)(2N+n+1)\\frac{(kT)^n}{n!} $$
                             where $$c$$ is the initiation rate, $$k$$ 
                             the elongation rate and $$T$$ is 
                             the residence time.

                             If you choose approximation_fit, the default the 
                             equation used is :
                             $$ G(x) = \\frac{T - x}{cT^2}
                             H(T - x), $$
                             where $$c$$ is the initiation rate, and $$T$$ is the residence time.

                             $$T=M/k$$ where $M$ is the RNA size (aa) and $k$ is the
                             elongation rate.
                            ''',
                             mathjax=True),
            ]),

            # Equation
            dbc.Col([
                html.H5("Equation"),

                dcc.Markdown('''
                You can change the equation. Symbols used are: 
                * x (fluorescence input), 
                * t (elongation rate) 
                * c (initiation rate). 
                
                If you need other symbols, please contact the admin.
                '''),
                dcc.Markdown(''' This part is not working (and it is not 
                useful) ''',
                             style={'color': 'red', }),
                dbc.Select(options=items_equation,
                            id="equation1",
                            className="mr-2",
                            style={"width": 350}, ),

                html.Br(),


                dbc.Input(id='equation2',
                           type='string',
                           value="",
                           placeholder="Enter equation",
                           style={'width': '350px'}),
                dbc.FormFeedback('Valid equation input!',
                                  type='valid'),
                dbc.FormFeedback('Invalid equation input!',
                                  type='invalid'),

                html.Br(),
                dbc.Button('Validate equation',
                            id='submit-button-equation2',
                            className="mr-2",
                            style={"width": "150px"}, ),


                # display the chosen equation and that it is valid
                dbc.Row([
                     html.Div(id='loading_equation_output2'),
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
        ]),
        # html.Br(),
        # Put here the chosen equation
        dbc.Row([
            dbc.Col(html.P("The chosen solver is: ", className="mb-0"),
                    width="auto"),
            dbc.Col(html.Div(id="choosen-solver1"), width="auto"),
        ]),
        # dbc.Row([
        #     dbc.Col(html.P("The chosen equation is :", className="mb-0"),
        #             width="auto"),
        #     dbc.Col(html.Div(id='equation_select0'), width="auto"),
        #
        # ]),

        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P("dt (sec)", style={"height": "auto",
                                              "margin-bottom": "auto"}),
                    dcc.Input(id='dt-param-vivo2', type='number', value=3),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(["Suntag length (aa)",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Length of the SunTag in amino "
                                "acid.",
                                target="faq_prot_length"),
                    dcc.Input(id='suntag-length-param-vivo2', type='number',
                              value=800),
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
                    dbc.Tooltip("Length of the protein in amino "
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
                    dcc.Input(id='missing_point_param_vivo2',
                              type='number',
                              value=5),
                ]),
            ]),
            dbc.Col([
                html.Div([
                    html.P(["Number of suntag",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Number of suntag repetition.",
                                target="faq_prot_length"),
                    dcc.Input(id='repetition-suntag-param-vivo2',
                              type='number',
                              value=32),
                ]),
            ]),

            dbc.Col([

            ]),

        ]),

        dbc.Row(
            [
                # dbc.Label(" "),
                dbc.Checklist(options=[{"label": "Force the analysis even if "
                                                 "time is not continuous.",
                                        "value": 0}],
                              id="switches_force_analysis2",
                              switch=False, ),
            ]),
        # Tick box to choose if we use the first dot for the analysis
        dbc.Row(
            [
                # dbc.Label(" "),
                dbc.Checklist(options=[{"label": "Include the first dot of "
                                                 "the autocorrelation "
                                                 "in the fit.",
                                        "value": 0}],
                              id="switches_first_dot2",
                              switch=False, ),
            ]),

        html.Br(),

        # Generate Button and Spinner Side by Side
        dbc.Row([
            # Analysis ONE track and ID of the track
            dbc.Col([
                html.H5("Display of one track",
                        style={"text-align": "center"}),
                dbc.Row([
                    html.Div([
                        html.P("id of the track to analyse",
                               style={"height": "auto",
                                      "margin-bottom": "auto"}),
                        dcc.Input(id='id_track2',
                                  type='number',
                                  value=0),
                    ]),
                ]),
                dbc.Row([
                    dbc.Col([
                        dbc.Button('Analyse and display ONE track',
                                   id='analyse_show_button2',
                                   className="mr-2",
                                   style={"width": "300px"}, ),
                    ], width="auto"),

                    dbc.Col([
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
            dbc.Col([
                html.H5("Analyse all tracks",
                        style={"text-align": "center"}),
                dbc.Row([
                    html.Div([
                        html.P("File name to save", style={"height": "auto",
                                                           "margin-bottom": "auto"}),
                        dcc.Input(id='save-results-name-vivo', type='text',
                                  value='datas_results'),
                    ]),
                ]),
                dbc.Row([
                    dbc.Col([
                        dbc.Button('Analyse ALL Tracks',
                                   id='start-analyze-btn-vivo',
                                   className="mr-2",
                                   style={"width": "300px"}, ),
                    ], width="auto"),

                    dbc.Col([
                        dbc.Spinner(
                            children=[html.Div(id="loading_analysis_vivo")],
                            size="sm", color="primary", type="border",
                            spinner_style={"margin-left": "10px"}
                        )
                    ], width="auto"),
                ], align="center", style={"margin-top": "10px",}),
                html.Div(id='analyze-output-vivo'),
                dcc.Download(id="download-csv"),
            ], style={'border': 'ridge',
                      'background-color': '#e6ffde'}),
        ]),

        html.Br(),
        html.Br(),

        dbc.Row([
            html.H5("Display of one track",
                    style={"text-align": "center"}),
            dcc.Graph(id='plot_results'),
        ]),

        html.Br(),
        html.Br(),

    )


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-vivo2', 'children'),
        Output("table-container2", "children"),
        Output('loading_data_vivo2', 'children'),
        Input('browse_directory_analyze_vivo2', 'contents'),
    )
    def browse_directory_analyze_vivo2(contents):
        df, output = upload_csv(contents, app, "csv_to_analyse")
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
        Output("choosen-solver1", "children"),
        # Output("equation_select0", "children", allow_duplicate=True),
        Input('choose-solver1', 'value'),
        prevent_initial_call=True,
    )
    def validate_solver2(value):
        app.data["solver"] = value
        # if value == "linear fit":
        #     app.data["equation_f"] = None
        #     app.data["equation_display"] = None
        return value #, str(app.data["equation_display"])


    @app.callback(
        Output("equation_select0", "children", allow_duplicate=True),
        Input('equation1', 'value'),
        prevent_initial_call=True,
    )
    def validate_input2(value):
        if value is None:
            return "None"
        else:
            v_bool, v_message, func_, expr_ = validate_equation(value)
            if value and v_bool:
                app.data["equation_f"] = func_
                app.data["equation_display"] = expr_

                return str(expr_)
            else:
                return ""


    @app.callback(
        Output('equation2', 'valid'),
        Output('equation2', 'invalid'),
        Output("loading_equation_output2", "children"),
        Output("equation_select0", "children"),
        Input('submit-button-equation2', 'n_clicks'),
        State('equation2', 'value'),
        prevent_initial_call=True,
    )
    def validate_input2(n_clicks, value):
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
        Output('plot_results', 'figure'),
        Output('loading_output_track', 'children'),
        Output('loading_output1', 'children'),
        Output('loading_output2', 'children'),
        Output('loading_output3', 'children'),
        Output('loading_track_plot', 'children'),
        Input('analyse_show_button2', 'n_clicks'),
        State('col_track2', 'value'), #0
        State('col_time2', 'value'), #1
        State('col_intensity2', 'value'), #2
        State('dt-param-vivo2', 'value'), #3
        State('prot-length-param-vivo2', 'value'), #4
        State('suntag-length-param-vivo2', 'value'), #5
        State('repetition-suntag-param-vivo2', 'value'), #6
        State('id_track2', 'value'), #7
        State("missing_point_param_vivo2", 'value'), #8
        State("switches_first_dot2", "value"), #9
        State("switches_force_analysis2", "value"), #10
    )
    def analyse_display_track(n_clicks, *params):
        figure = make_subplots(rows=3,
                               cols=1,
                               subplot_titles=["track profile",
                                               "autocorrelation",
                                               "residuals"
                                               ])
        if n_clicks:
            try:
                str_output_force = ""
                # Read csv file
                df = app.data['csv_to_analyse']

                df.rename(columns={params[0]: 'TRACK_ID',
                                      params[1]: 'FRAME',
                                      params[2]: 'MEAN_INTENSITY_CH1',
                                      },
                             inplace=True)

                if int(params[7]) not in np.unique(df["TRACK_ID"]):
                    return figure, "", "This ID does not exist.", "", "", None

                dt = float(params[3])
                prot_length = int(params[4])
                suntag_length = int(params[5])
                repetition_suntag = int(params[6])
                datas2 = df[(df["TRACK_ID"] == int(params[7]))]
                first_dot = bool(params[-2])
                force_analysis = bool(params[-1])

                # Check solver/equation parameters are present
                # if valid, track can be analysed
                if app.data["solver"] == "exact_fit":
                    method = "exact"

                elif app.data["solver"] == "approximation fit":
                    method = "approx"
                    # TODO : check equation are selected

                elif app.data["solver"] == "linear fit":
                    method = "linear"

                valid, x, y, x_fix, y_fix = check_track_validity(datas2,
                                            int(params[7]),
                                            normalise_intensity=1,
                                            delta_t = dt,
                                            rtol=1e-1,
                                            nb_missing_point=int(params[-3]),
                                            )



                if valid:
                    (x_auto,
                     y_auto,
                     k, c,
                     elongation_r,
                     translation_init_r,
                     perr) = single_track_analysis(x_fix,
                                                   y_fix,
                                                   delta_t=dt,
                                                   protein_size=prot_length,
                                                   suntag_size=suntag_length,
                                                   repetition_suntag=repetition_suntag,
                                                   mm=None,
                                                   normalise_auto=True,
                                                   method=method,
                                                   first_dot=first_dot,
                                                   simulation=False,
                                                   func_=app.data["equation_f"]
                                                   )
                else:
                    if force_analysis:
                        (x_auto,
                         y_auto,
                         k,c,
                         elongation_r,
                         translation_init_r,
                         perr) = single_track_analysis(x_fix,
                                                       y_fix,
                                                       delta_t=dt,
                                                       protein_size=prot_length,
                                                       mm=None,
                                                       normalise_auto=True,
                                                       method=method,
                                                       first_dot=first_dot,
                                                       simulation=False,
                                                       func_=app.data[
                                                           "equation_f"]
                                                       )
                        str_output_force = f"Analysis has been forced!"
                    else:
                        return (figure, "This track can't be analysed, "
                                        "there is too much missing point. "
                                        "If you want to force the analysis, "
                                        "you can either increase the missing "
                                        "point value or tick the force "
                                        "analysis box.",
                                "", "","", None)




                # plot track profile
                figure.add_trace(go.Scatter(x=x*dt,
                                            y=y,
                                            mode="lines",
                                            name="Signal",
                                            line_color="#1A51DB"), #Blue
                                            row=1,
                                            col=1,
                                            )

                figure.add_trace(go.Scatter(x=x_fix*dt,
                                            y=y_fix,
                                            mode="lines+markers",
                                            name="Signal correct",
                                            line_color="#DBA41A"), #Orange
                                 row=1,
                                 col=1,
                                 )
                figure.update_xaxes(title_text='Time (sec)', row=1, col=1)
                figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

                # plot track profile
                figure.add_trace(go.Scatter(x=x_auto,
                                            y=y_auto,
                                            mode="lines",
                                            name="Autocorrelation",
                                            line_color="#000000"), #Black
                                            row=2,
                                            col=1),

                figure.add_trace(go.Scatter(x=x_auto,
                                            y=y_auto,
                                            mode="markers",
                                            name="Autocorrelation",
                                            line_color="#000000"),  # Black
                                 row=2,
                                 col=1),
                if method == "exact":
                    figure.add_trace(
                        go.Scatter(x=x_auto[:int(len(x_auto) / 2)],
                                   y=fit_function_exact(x_auto,
                                                        k, c, repetition_suntag)[
                                       :int(len(x_auto) / 2)],
                                   mode="lines",
                                   name="Fit",
                                   line_color="#B80909"),  # Red
                        row=2,
                        col=1),

                    figure.add_trace(
                        go.Scatter(x=x_auto[:int(len(x_auto) / 2)],
                                   y=fit_function_exact(x_auto,
                                                        k, c, repetition_suntag)[
                                       :int(len(x_auto) / 2)],
                                   mode="markers",
                                   name="Fit",
                                   line_color="#B80909"),  # Red
                        row=2,
                        col=1),

                elif method =="approx" :
                    figure.add_trace(go.Scatter(x=x_auto[:int(len(x_auto)/2)],
                                                y=fit_function(x_auto,
                                                               k, c)[:int(len(
                                                    x_auto)/2)],
                                                mode="lines",
                                                name="Fit",
                                                line_color="#B80909"), #Red
                                     row=2,
                                     col=1),

                    figure.add_trace(go.Scatter(x=x_auto[:int(len(x_auto) / 2)],
                                                y=fit_function(x_auto,
                                                               k, c)[:int(len(
                                                    x_auto) / 2)],
                                                mode="markers",
                                                name="Fit",
                                                line_color="#B80909"),  # Red
                                     row=2,
                                     col=1),

                else :
                    figure.add_trace(go.Scatter(x=x_auto[:int(-c/k)+1],
                                                y=(x_auto*k+c)[:int(-c/k)+1],
                                                mode="lines",
                                                name="Fit",
                                                line_color="#B80909"), #Red
                                     row=2,
                                     col=1),

                    figure.add_trace(go.Scatter(x=x_auto[:int(-c/k)+1],
                                                y=(x_auto*k+c)[:int(-c/k)+1],
                                                mode="markers",
                                                name="Fit",
                                                line_color="#B80909"),  # Red
                                     row=2,
                                     col=1),


                figure.update_xaxes(title_text='Time delay (tau)',
                                    row=2,
                                    col=1)
                figure.update_yaxes(title_text='G(tau)', row=2, col=1)
                figure.update_layout(width=1000, height=800, )

                # plot residuals
                figure.add_trace(go.Scatter(x=x_auto[:int(len(x_auto) / 2)],
                                            y=np.repeat(0, len(x_auto[:int(
                                                len(x_auto) / 2)])),
                                            mode="lines",
                                            name="",
                                            line_color="#aeb6c2"),  # grey
                                 row=3,
                                 col=1)

                figure.add_trace(go.Scatter(x=x_auto,
                                            y=y_auto[:int(len(x_auto)/2)]-fit_function(x_auto,
                                                           prot_length/elongation_r, 1/translation_init_r)[:int(len(x_auto)/2)],
                                            mode="markers",
                                            name="residuals",
                                            line_color = "#000000"),  # black
                                 row=3,
                                 col=1)



                figure.update_xaxes(title_text='X data',
                                    row=3,
                                    col=1)
                figure.update_yaxes(title_text='Delta', row=3, col=1)
                figure.update_layout(width=1000, height=800, )

                str_output1 = (f"elongation rate: "
                              f"{elongation_r:.2f} aa/sec ")
                str_output2 = (f"initiation rate: "
                              f"{translation_init_r:.2f} rib/sec")

                str_output3 = (f"perr: "
                               f"{perr} ")
                return figure, str_output_force, str_output1, str_output2, str_output3, None
            except Exception as e:
                print(e)
                return {
                    'data': [],
                    'layout': go.Layout(title='Error', xaxis={'title': 'Time'},
                                        yaxis={'title': 'Fluorescence'})
                },"",  str(e),"", "", None
        return figure, "",  "", "", "", None


    @app.callback(
        Output('analyze-output-vivo', 'children'),
        Output('loading_analysis_vivo', 'children'),
        Output('download-csv', 'data'),
        Input('start-analyze-btn-vivo', 'n_clicks'),

        State('col_track2', 'value'), #0
        State('col_time2', 'value'), #1
        State('col_intensity2', 'value'), #2
        State('dt-param-vivo2', 'value'), #3
        State('prot-length-param-vivo2', 'value'), #4
        State('suntag-length-param-vivo2', 'value'), #5
        State('repetition-suntag-param-vivo2', 'value'), #6
        State("missing_point_param_vivo2", 'value'), #7
        State("switches_first_dot2", "value"), #8
        State("switches_force_analysis2", "value"), #9
        State('save-results-name-vivo', 'value'), #10
        # State('checkbox_simu', 'value') #9
    )
    def start_analyze_all_tracks(n_clicks, *params):

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
                suntag_length = int(params[5])
                repetition_suntag = int(params[6])
                nb_missing_point = int(params[7])

                first_dot = bool(params[8])
                force_analysis = bool(params[9])
                # check_simu = 'checked' in params[9]

                # Check solver/equation parameters are present
                # if valid, track can be analysed
                if app.data["solver"] == "exact_fit":
                    method = "exact"
                elif app.data["solver"] == "approximation fit":
                    method = "approx"
                    # TODO : check equation are selected

                elif app.data["solver"] == "linear fit":
                    method = "linear"

                ids_track = np.unique(df["TRACK_ID"])
                first_time = True
                # Analyse all tracks and save it
                for i in ids_track:
                    print(i)
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
                         k,c,
                         elongation_r,
                         translation_init_r,
                         perr) = single_track_analysis(x_fix,
                                                       y_fix,
                                                       delta_t=dt,
                                                       protein_size=prot_length,
                                                       suntag_size=suntag_length,
                                                       repetition_suntag=repetition_suntag,
                                                       mm=None,
                                                       normalise_auto=True,
                                                       method=method,
                                                       first_dot=first_dot,
                                                       simulation=False,
                                                       func_=app.data["equation_f"]
                                                       )
                        print(k, c, elongation_r, translation_init_r)
                        if first_time:
                            results = pd.DataFrame({
                                                    "id": i,
                                                    "dt": dt,
                                                    "length": len(x_fix),
                                                    "k":k,
                                                    "c":c,
                                                    "elongation_r": elongation_r,
                                                    "init_translation_r": translation_init_r,
                                                    "perr0":perr[0],
                                                    "perr1": perr[1],
                                                    "comment": comment},
                                                   index=[0])
                            first_time = False

                        else:
                            results = pd.concat([results,
                                                 pd.DataFrame(
                                                     {
                                                        "id": i,
                                                        "dt": dt,
                                                        "length": len(x_fix),
                                                        "k": k,
                                                        "c": c,
                                                        "elongation_r": elongation_r,
                                                        "init_translation_r": translation_init_r,
                                                        "perr0":perr[0],
                                                        "perr1": perr[1],
                                                        "comment":comment},
                                                     index=[0])
                                                 ], ignore_index=True)
                    else:
                        if first_time:
                            results = pd.DataFrame(
                                {
                                 "id": i,
                                 "dt": np.nan,
                                 "length": len(x_fix),
                                 "k": np.nan,
                                 "c": np.nan,
                                 "elongation_r": np.nan,
                                 "init_translation_r": np.nan,
                                 "perr0": np.nan,
                                 "perr1": np.nan,
                                 "comment":"cant be analysed"},
                                index=[0])
                            first_time = False

                        else:
                            results = pd.concat([results,
                                                 pd.DataFrame(
                                                     {
                                                         "id": i,
                                                         "dt": np.nan,
                                                         "length": len(x_fix),
                                                         "k": np.nan,
                                                         "c": np.nan,
                                                         "elongation_r":
                                                             np.nan,
                                                         "init_translation_r": np.nan,
                                                         "perr0": np.nan,
                                                         "perr1": np.nan,
                                                     "comment":"cant be "
                                                               "analysed"},
                                                     index=[0])
                                                 ], ignore_index=True)



                output_path = params[-1] + ".csv"
                results.to_csv(output_path, index=False)

                return "Analysis completed and saved successfully!", None, dcc.send_file(output_path)
            except Exception as e:
                return f"Error: {str(e)}", None, None
        raise PreventUpdate
