import numpy as np
import pandas as pd

from dash import dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from ..tabs.app_function import (upload_csv)
from ..analysis.analysis_track import (single_track_analysis,
                                       check_track_validity)
from ..analysis.fit_funtions import (fit_function_exact,
                                     fit_function_approx,
                                     fit_function_epitope)
from ..plots.plots import fig_analyse_track


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-vivo2', 'children'),
        Output("table-container2", "children"),
        Output('loading_data_vivo2', 'children'),
        Input('browse_directory_analyze_vivo2', 'contents'),
        prevent_initial_call=True,
    )
    def browse_directory_analyze_vivo2(contents):
        df, output = upload_csv(contents, app, "csv_to_analyse")
        if df is None:
            return output, None, None

        table = dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{"name": i, "id": i} for i in
                     df.columns],
            page_size=10,
            style_table={'width': '800px', 'overflowX': 'auto'},
        )
        return None, table, None

    @app.callback(
        Output("choosen-solver1", "children"),
        Input('choose-solver1', 'value'),
        prevent_initial_call=True,
    )
    def validate_solver2(value):
        app.data["solver"] = value
        return value


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
        State("switches_force_analysis2", "value"), #9
    )
    def analyse_display_track(n_clicks, *params):

        if n_clicks:
            try:
                str_output_force = ""
                # get table and rename columns if needed
                df = app.data['csv_to_analyse']
                df.rename(columns={params[0]: 'TRACK_ID',
                                   params[1]: 'FRAME',
                                   params[2]: 'MEAN_INTENSITY_CH1',
                                   },
                             inplace=True)

                if int(params[7]) not in np.unique(df["TRACK_ID"]):
                    return go.Figure(), "", "This ID does not exist.", "", "", None

                dt = float(params[3])
                prot_length = int(params[4])
                suntag_length = int(params[5])
                repetition_suntag = int(params[6])
                datas2 = df[(df["TRACK_ID"] == int(params[7]))]
                force_analysis = bool(params[-1])

                # Check solver
                if app.data["solver"] == "Exact equation":
                    method = "exact"
                elif app.data["solver"] == "Approximate equation":
                    method = "approx"
                elif app.data["solver"] == "Approximate epitope":
                    method = "epitope"


                valid, x, y, x_fix, y_fix = check_track_validity(datas2,
                                            int(params[7]),
                                            normalise_intensity=1,
                                            delta_t = dt,
                                            rtol=1e-1,
                                            nb_missing_point=int(params[-3]),
                                            )

                if valid or force_analysis:
                    if force_analysis:
                        str_output_force = f"Analysis has been forced!"
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
                                                   simulation=False,
                                                   )
                else:
                    return (go.Figure(), "This track can't be analysed, "
                                    "there is too much missing point. "
                                    "If you want to force the analysis, "
                                    "you can either increase the missing "
                                    "point value or tick the force "
                                    "analysis box.",
                            "", "","", None)

                N = repetition_suntag
                M = prot_length/(suntag_length/repetition_suntag)
                if method == "exact":
                    y_fit = fit_function_exact(x_auto, k, c, N, M)
                elif method == "approx":
                    y_fit = fit_function_approx(x_auto, k, c)
                elif method == "epitope":
                    y_fit = fit_function_epitope(x_auto, k, c, N)

                figure = fig_analyse_track(x, y,
                                           x_fix, y_fix,
                                           x_auto, y_auto,
                                           y_fit,
                                           dt)


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
        return go.Figure(), "",  "", "", "", None


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
        State("switches_force_analysis2", "value"), #8
        State('save-results-name-vivo', 'value'), #9
        # State('checkbox_simu', 'value') #9
    )
    def start_analyze_all_tracks(n_clicks, *params):

        if n_clicks:
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
                force_analysis = bool(params[8])
                # check_simu = 'checked' in params[9]

                # Check solver
                if app.data["solver"] == "Exact equation":
                    method = "exact"
                elif app.data["solver"] == "Approximate equation":
                    method = "approx"
                elif app.data["solver"] == "Approximate epitope":
                    method = "epitope"

                ids_track = np.unique(df["TRACK_ID"])
                first_time = True
                # Analyse all tracks and save it
                for i in ids_track:
                    print(i)

                    k = np.nan
                    c = np.nan
                    elongation_r = np.nan
                    translation_init_r = np.nan
                    perr0 = np.nan
                    perr1 = np.nan
                    comment = ""

                    datas2 = df[(df["TRACK_ID"] == i)]

                    (valid, x, y, x_fix, y_fix) = check_track_validity(datas2,
                                                   i,
                                                   normalise_intensity=1,
                                                   delta_t=dt,
                                                   rtol=1e-1,
                                                   nb_missing_point=nb_missing_point,
                                                   )
                    length = len(x_fix)

                    if valid or force_analysis:
                        if force_analysis:
                            comment = "analysis forced"
                        (x_auto,
                         y_auto,
                         k,c,
                         elongation_r,
                         translation_init_r,
                         [perr0, perr1]) = single_track_analysis(x_fix,
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

                    # Populate the dataframe
                    if first_time:
                        results = pd.DataFrame({
                                                "id": i,
                                                "dt": dt,
                                                "length": length,
                                                "k":k,
                                                "c":c,
                                                "elongation_r": elongation_r,
                                                "init_translation_r": translation_init_r,
                                                "perr0":perr0,
                                                "perr1": perr1,
                                                "comment": comment},
                                               index=[0])
                        first_time = False

                    else:
                        results = pd.concat([results,
                                             pd.DataFrame(
                                                 {
                                                    "id": i,
                                                    "dt": dt,
                                                    "length": length,
                                                    "k": k,
                                                    "c": c,
                                                    "elongation_r": elongation_r,
                                                    "init_translation_r": translation_init_r,
                                                    "perr0":perr0,
                                                    "perr1": perr1,
                                                    "comment":comment},
                                                 index=[0])
                                             ], ignore_index=True)

                output_path = params[-1] + ".csv"
                results.to_csv(output_path, index=False)

                return "Analysis completed and saved successfully!", None, dcc.send_file(output_path)
            except Exception as e:
                return f"Error: {str(e)}", None, None
        raise PreventUpdate
