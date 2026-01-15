import numpy as np

from dash import  html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate

from .utils import generate_table, generate_table_selectable
from ..tabs.app_function import (upload_csv)

from ..analysis.movement import msd_calculation, rsd_calculation

## Callbacks
def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-msd', 'children'),
        Output("table_container_msd", "children"),
        Output('loading_data_msd', 'children'),
        Input('browse_directory_msd', 'contents'),
    )
    def browse_directory_msd(contents):
        df, output = upload_csv(contents, app, "csv_msd")
        if df is None:
            return output, None, None

        table = generate_table(df, 10)

        return None, table, None

    @app.callback(
        Output('output_MSD', 'children'),
        Output('download_csv_msd', 'data'),
        Input('btn_calculate_MSD', 'n_clicks'),
        State('col_ID', 'value'),
        State('col_x', 'value'),
        State('col_y', 'value'),
        State('col_z', 'value'),
    )
    def calculate_msd(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """

        if n_clicks:
            if "csv_msd" not in app.data:
                return "You need to upload a file first", None

            df = app.data['csv_msd']
            df.rename(columns={params[0]: 'TRACK_ID',
                               params[1]: 'POSITION_X',
                               params[2]: 'POSITION_Y',
                               params[3]: 'POSITION_Z',
                               },
                      inplace=True)

            try:

                # MSD
                msd_val = []
                for i in np.unique(df["TRACK_ID"]):
                    x = df[df["TRACK_ID"] == i]["POSITION_X"]
                    y = df[df["TRACK_ID"] == i]["POSITION_Y"]
                    z = df[df["TRACK_ID"] == i]["POSITION_Z"]
                    msd = msd_calculation(x,y,z)
                    msd_val = np.concatenate([msd_val, msd])
                df["msd"] = msd_val

                # RSD
                rsd_val = []
                for i in np.unique(df["TRACK_ID"]):
                    x = df[df["TRACK_ID"] == i]["POSITION_X"].to_numpy()
                    y = df[df["TRACK_ID"] == i]["POSITION_Y"].to_numpy()
                    z = df[df["TRACK_ID"] == i]["POSITION_Z"].to_numpy()
                    rsd = rsd_calculation(x, y, z)
                    rsd_val = np.concatenate([rsd_val, rsd])
                df["rsd"] = rsd_val

                output_path = "result.csv"
                df.to_csv(output_path, index=False)
                return "", dcc.send_file(output_path)
            except Exception as e:
                print(e)
                return "Problem", None
        raise PreventUpdate
