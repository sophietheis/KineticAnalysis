import logging

from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate

from .utils import generate_table_selectable, empty_error_figure
from ..session_store import (get_session_data,
                             has_session_data,
                             set_session_data)
from ..plots.plots import plot_ribosome, plot_ribosome_density
from ..tabs.app_function import (upload_csv)

from ..analysis.analyse_density import calculate_ribosome_density

logger = logging.getLogger(__name__)


def register_callbacks(app):
    @app.callback(
        Output('selected-file-output-single_prot', 'children'),
        Output("table_container_single", "children"),
        Output('loading_data_single_prot', 'children'),
        Input('browse_directory_single_prot', 'contents'),
        State('session-id', 'data'),
    )
    def browse_directory_single_prot(contents, session_id):
        df, output = upload_csv(contents, session_id, "csv_fluo_single")
        if df is None:
            return output, None, None

        table = generate_table_selectable(df, max_rows=10, width="300px",
                                          **{"id": "table_single_prot"})

        return None, table, None

    @app.callback(
        Output('table_single_prot', 'style_data_conditional'),
        Input('table_single_prot', 'selected_columns'),
        State('session-id', 'data'),
    )
    def single_prot_select_name(selected_columns, session_id):
        if len(selected_columns) == 0:
            return None

        set_session_data(session_id,
                         "single_prot_column_intensity",
                         selected_columns[0])
        return [{'if': {'column_id': i},
                 'background_color': '#D2F3FF'
                 } for i in selected_columns]

    @app.callback(
        Output('selected-file-output-polysome', 'children'),
        Output("table_container_polysome", "children"),
        Output('loading_data_polysome', 'children'),
        Input('browse_directory_polysome', 'contents'),
        State('session-id', 'data'),
    )
    def browse_directory_polysome(contents, session_id):
        df, output = upload_csv(contents, session_id, "csv_fluo_polysome")
        if df is None:
            return output, None, None

        table = generate_table_selectable(df, max_rows=10, width="300px",
                                          **{"id": "table_polysome", })
        return None, table, None

    @app.callback(
        Output('table_polysome', 'style_data_conditional'),
        Input('table_polysome', 'selected_columns'),
        State('session-id', 'data'),
    )
    def polysome_select_name(selected_columns, session_id):
        if len(selected_columns) == 0:
            return None

        set_session_data(session_id,
                         "polysome_column_intensity",
                         selected_columns[0])
        return [{'if': {'column_id': i},
                 'background_color': '#D2F3FF'
                 } for i in selected_columns]

    @app.callback(
        Output('ribosome-plot', 'figure'),
        Output('ribosome-density-plot', 'figure'),
        Output('output_ribosome_density', 'children'),
        Output('download_csv', 'data'),
        Input('btn_calculate_ribosome_density', 'n_clicks'),
        State('session-id', 'data'),
        State('param_prot_length_rib', 'value'),  #0
        State('param_suntag_length_rib', 'value'),  #1
    )
    def calculate_density(n_clicks, session_id, *params):
        """
        This function generate and plot an example for the simulation.
        """

        if n_clicks:
            if not has_session_data(session_id,
                                    "polysome_column_intensity"):
                return (None,
                        None,
                        "Please select a column in polysome dataframe",
                        None)
            if not has_session_data(session_id,
                                    "single_prot_column_intensity"):
                return (None,
                        None,
                        "Please select a column in single prot dataframe",
                        None)

            try:
                L_poi = float(params[0])
                L_tag = float(params[1])

                m_intensity_single, result = calculate_ribosome_density(
                    get_session_data(session_id, "csv_fluo_single"),
                    get_session_data(session_id, "csv_fluo_polysome"),
                    get_session_data(session_id, "single_prot_column_intensity"),
                    get_session_data(session_id, "polysome_column_intensity"),
                    L_poi,
                    L_tag)

                mean_single = get_session_data(session_id, "csv_fluo_single")[get_session_data(session_id, "single_prot_column_intensity")].mean()
                std_single = get_session_data(session_id, "csv_fluo_single")[get_session_data(session_id, "single_prot_column_intensity")].std()
                mean_polysome = result["INTENSITY"].mean()
                std_polysome = result["INTENSITY"].std()
                mean_rib = result["ribosome_density"].mean()
                std_rib = result["ribosome_density"].std()
                output_string = html.P([
                    f"Mean single protein intensity : {mean_single:.2f}",
                    html.Br(),
                    f"STD single protein intensity : {std_single:.2f}",
                    html.Br(),
                    f"Mean polysome intensity : {mean_polysome:.2f}",
                    html.Br(),
                    f"STD polysome intensity : {std_polysome:.2f}",
                    html.Br(),
                    f"Mean ribosome density : {mean_rib:.2f} rib/aa, {1/mean_rib:.2f} aa/rib,",
                    html.Br(),
                    f"STD ribosome density : {std_rib:.2f} rib/aa, {1/std_rib:.2f} aa/rib"
                    ])

                output_path = "result.csv"
                result.to_csv(output_path, index=False)
                figure = plot_ribosome(get_session_data(session_id, "csv_fluo_single")[get_session_data(session_id, "single_prot_column_intensity")],
                                       get_session_data(session_id, "csv_fluo_polysome")[get_session_data(session_id, "polysome_column_intensity")],
                                       result
                                       )

                figure_density = plot_ribosome_density(L_poi, L_tag, mean_rib)
                return (figure,
                        figure_density,
                        output_string,
                        dcc.send_file(output_path))
            except Exception as e:
                logger.exception("Failed to calculate ribosome density", e)
                return (empty_error_figure(),
                        empty_error_figure(),
                        "Problem",
                        None)
        raise PreventUpdate
