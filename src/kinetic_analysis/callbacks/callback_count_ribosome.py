import numpy as np

from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate

from .utils import generate_table_selectable, empty_error_figure
from ..plots.plots import plot_ribosome
from ..tabs.app_function import (upload_csv)

from ..analysis.analyse_density import calculate_ribosome_density


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-single_prot', 'children'),
        Output("table_container_single", "children"),
        Output('loading_data_single_prot', 'children'),
        Input('browse_directory_single_prot', 'contents'),
    )
    def browse_directory_single_prot(contents):
        df, output = upload_csv(contents, app, "csv_fluo_single")
        if df is None:
            return output, None, None

        table = generate_table_selectable(df, max_rows=10, width="300px",
                                          **{"id": "table_single_prot"})

        return None, table, None

    @app.callback(
        Output('table_single_prot', 'style_data_conditional'),
        Input('table_single_prot', 'selected_columns'),
    )
    def single_prot_select_name(selected_columns):
        if len(selected_columns) == 0:
            return None

        app.data["single_prot_column_intensity"] = selected_columns[0]
        return [{'if': {'column_id': i},
                 'background_color': '#D2F3FF'
                 } for i in selected_columns]

    @app.callback(
        Output('selected-file-output-polysome', 'children'),
        Output("table_container_polysome", "children"),
        Output('loading_data_polysome', 'children'),
        Input('browse_directory_polysome', 'contents'),
    )
    def browse_directory_polysome(contents):
        df, output = upload_csv(contents, app, "csv_fluo_polysome")
        if df is None:
            return output, None, None

        table = generate_table_selectable(df, max_rows=10, width="300px",
                                          **{"id": "table_polysome", })
        return None, table, None

    @app.callback(
        Output('table_polysome', 'style_data_conditional'),
        Input('table_polysome', 'selected_columns'),
    )
    def polysome_select_name(selected_columns):
        if len(selected_columns) == 0:
            return None
        app.data["polysome_column_intensity"] = selected_columns[0]
        return [{'if': {'column_id': i},
                 'background_color': '#D2F3FF'
                 } for i in selected_columns]

    @app.callback(
        Output('ribosome-plot', 'figure'),
        Output('output_ribosome_density', 'children'),
        Output('download_csv', 'data'),
        Input('btn_calculate_ribosome_density', 'n_clicks'),
        State('param_prot_length_rib', 'value'),  #0
        State('param_suntag_length_rib', 'value'),  #1
    )
    def calculate_density(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """

        if n_clicks:
            if "polysome_column_intensity" not in app.data.keys():
                return None, "Please select a column in polysome dataframe", None
            if "single_prot_column_intensity" not in app.data.keys():
                return None, "Please select a column in single prot dataframe", None
            try:
                L_poi = float(params[0])
                L_tag = float(params[1])

                m_intensity_single, result = calculate_ribosome_density(
                    app.data["csv_fluo_single"],
                    app.data["csv_fluo_polysome"],
                    app.data["single_prot_column_intensity"],
                    app.data["polysome_column_intensity"],
                    L_poi,
                    L_tag)

                output_string = html.P([
                    # f"Mean single protein intensity : {np.round(app.data["csv_fluo_single"][app.data["single_prot_column_intensity"]].mean(), 2)}",
                    # f"STD single protein intensity : {np.round(app.data["csv_fluo_single"][app.data["single_prot_column_intensity"]].std(), 2)}",
                    f"Mean single protein intensity : {np.round(m_intensity_single, 2)}",
                    html.Br(),
                    # f"Mean polysome intensity : {np.round(result["INTENSITY"].mean(), 2)}",
                    # f"STD polysome intensity : {np.round(result["INTENSITY"].std(), 2)}",
                    html.Br(),
                    # f"Mean ribosome density : {np.round(result["ribosome_density"].mean(), 2)} rib/aa"
                    # f"STD ribosome density : {np.round(result["ribosome_density"].std(), 2)} rib/aa"
                    ])

                output_path = "result.csv"
                result.to_csv(output_path, index=False)
                figure = plot_ribosome(app.data["csv_fluo_single"][app.data["single_prot_column_intensity"]],
                                       app.data["csv_fluo_polysome"][app.data["polysome_column_intensity"]],
                                       L_poi, 
                                       L_tag,)

                return figure, output_string, dcc.send_file(output_path)
            except Exception as e:
                print(e)
                return empty_error_figure(), "Problem", None
        raise PreventUpdate
