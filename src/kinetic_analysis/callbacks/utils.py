from dash import dash_table

from dash import Input, Output, State

import pandas as pd

import plotly.graph_objs as go

from ..session_store import get_session_data
from ..tabs.app_function import upload_csv


def generate_table(dataframe, max_rows=10, width="800px", **kwargs):
    table = dash_table.DataTable(
        data=dataframe.to_dict('records'),
        columns=[{"name": i, "id": i} for i in
                 dataframe.columns],
        page_size=max_rows,
        style_table={'width': width, 'overflowX': 'auto'},
        **kwargs,
    )
    return table


def generate_table_selectable(dataframe, max_rows=10, width="800px", **kwargs):
    table = dash_table.DataTable(
        data=dataframe.to_dict('records'),
        columns=[{"name": i, "id": i, "selectable": True} for i in
                 dataframe.columns],
        page_size=max_rows,
        column_selectable="single",
        selected_columns=[],
        style_table={'width': width, 'overflowX': 'auto'},
        **kwargs,
    )
    return table


def _empty_selectable_table(table_id):
    """
    """
    return generate_table_selectable(pd.DataFrame(),
                                     max_rows=10,
                                     width="300px",
                                     id=table_id)


def resolve_solver_method(session_id):
    solver = get_session_data(session_id, "solver", default="Exact equation")
    return {
        "Epitope equation": "epitope",
        "Exact equation": "exact",
        "Post-epitope equation": "approx",
    }.get(solver,  "exact")


def empty_error_figure():
    return {
        "data": [],
        "layout": go.Layout(
            title='Error',
            xaxis={'visible': False},
            yaxis={'visible': False},
            annotations=[
                {
                    'text': "Error: No data to display",
                    'xref': 'paper',
                    'yref': 'paper',
                    'showarrow': False,
                    'font': {
                        'size': 20
                    }
                }
            ]
        )
    }


def register_upload_callback(app, *, contents_id, output_id, table_container_id,
                             loading_id, session_data_key, table_builder,
                             extra_outputs=None, prevent_initial_call=None):
    extra_outputs = extra_outputs or []
    outputs = [
        Output(output_id, 'children'),
        Output(table_container_id, 'children'),
        Output(loading_id, 'children'),
    ] + [out for out, _ in extra_outputs]

    @app.callback(
        *outputs,
        Input(contents_id, 'contents'),
        State('session-id', 'data'),
        prevent_initial_call=prevent_initial_call,
    )
    def _upload_callback(contents, session_id):
        df, output = upload_csv(contents, session_id, session_data_key)
        if df is None:
            return (output, None, None) + tuple(None for _ in extra_outputs)

        table = table_builder(df)
        extra_values = tuple(compute(df) for _, compute in extra_outputs)
        return (None, table, None) + extra_values
