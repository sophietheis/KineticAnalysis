import numpy as np

from dash import dcc, Input, Output, State
from dash.exceptions import PreventUpdate

from .utils import generate_table
from ..tabs.app_function import (upload_csv)
from ..analysis.combination import combine_tracks_df, combine_tracks_df_ensemble


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-combine', 'children'),
        Output("table-container_combine", "children"),
        Output('loading_data_combine', 'children'),
        Output('nb_tracks_init', 'children'),
        Input('browse_directory_combine', 'contents'),
        prevent_initial_call=True,
    )
    def browse_directory_combine(contents):
        df, output = upload_csv(contents, app, "csv_to_combine")
        if df is None:
            return output, None, None, None

        table = generate_table(df, 10)
        nb_track = len(np.unique(df["TRACK_ID"]))

        return None, table, None, nb_track

    @app.callback(
        Output('combine-tracks-output', 'children'),
        Output('download-csv3', 'data'),
        Output("loading_combination", "children"),
        Output("complete2", "data"),
        Input("combine-tracks-btn", "n_clicks"),
        Input("start2", "data"),
        State('nb_tracks', "value"),
        State('nb_new_tracks', "value"),

    )
    def create_tracks(n_clicks, data, *params):
        if n_clicks:
            try:
                df = app.data['csv_to_combine']

                new_df = combine_tracks_df(df, int(params[1]), int(params[0]))

                output_path = "new_combined_tracks.csv"
                new_df.to_csv(output_path, index=False)
                return ("New combined tracks generated successfully !!",
                        dcc.send_file(output_path),
                        "Generate tracks",
                        data)
            except Exception as e:
                return (f"Error: {str(e)}",
                        None,
                        "Generate tracks",
                        data)
        raise PreventUpdate


    @app.callback(
        Output('combine-tracks-gtau-output', 'children'),
        Output('download-csv4', 'data'),
        Output("loading_combination_gtau", "children"),
        Output("complete3", "data"),
        Input("combine-tracks-gtau-btn", "n_clicks"),
        Input("start3", "data"),
        State('nb_tracks', "value"),
        State('nb_new_tracks', "value"),
    )
    def create_tracks_gtau(n_clicks, data, *params):
        if n_clicks:
            try:
                df = app.data['csv_to_combine']

                new_df = combine_tracks_df_ensemble(df, int(params[1]), int(params[0]))

                output_path = "new_combined_gtau_tracks.csv"
                new_df.to_csv(output_path, index=False)
                return ("New combined tracks generated successfully !!",
                        dcc.send_file(output_path),
                        "Generate tracks",
                        data)
            except Exception as e:
                return (f"Error: {str(e)}",
                        None,
                        "Generate tracks",
                        data)
        raise PreventUpdate