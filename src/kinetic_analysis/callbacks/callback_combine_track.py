import logging
import numpy as np

from dash import dcc, Input, Output, State
from dash.exceptions import PreventUpdate

from .utils import generate_table
from ..session_store import get_session_data
from ..tabs.app_function import (upload_csv)
from ..analysis.combination import combine_tracks_df

logger = logging.getLogger(__name__)


def register_callbacks(app):
    @app.callback(
        Output('combine_output_file', 'children'),
        Output("combine_table_container", "children"),
        Output('combine_loading_data', 'children'),
        Output('nb_tracks_init', 'children'),
        Input('combine_browse_directoy', 'contents'),
        State('session-id', 'data'),
        prevent_initial_call=True,
    )
    def combine_browse_directoy(contents, session_id):
        df, output = upload_csv(contents, session_id, "csv_to_combine")
        if df is None:
            return output, None, None, None

        table = generate_table(df, 10)
        nb_track = len(np.unique(df["TRACK_ID"]))

        return None, table, None, nb_track

    @app.callback(
        Output('combine_output_tracks', 'children'),
        Output('combine_download_csv', 'data'),
        Output("combine_loading", "children"),
        Output("combine_complete", "data"),
        Input("combine_button_combine", "n_clicks"),
        Input("combine_start", "data"),
        State('session-id', 'data'),
        State('combine_nb_tracks', "value"),
        State('combine_nb_new_tracks', "value"),

    )
    def create_tracks(n_clicks, data, session_id, *params):
        if n_clicks:
            try:
                df = get_session_data(session_id, 'csv_to_combine')

                new_df = combine_tracks_df(df, int(params[1]), int(params[0]))

                output_path = "new_combined_tracks.csv"
                new_df.to_csv(output_path, index=False)
                return ("New combined tracks generated successfully !!",
                        dcc.send_file(output_path),
                        "Generate tracks",
                        data)
            except Exception as e:
                logger.exception("Failed to combine tracks")
                return (f"Error: {str(e)}",
                        None,
                        "Generate tracks",
                        data)
        raise PreventUpdate
