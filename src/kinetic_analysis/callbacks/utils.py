from dash import dash_table

import plotly.graph_objs as go

from ..session_store import get_session_data


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

def resolve_solver_method(session_id):
    solver = get_session_data(session_id, "solver", default="Exact equation")
    return {
        "Exact equation": "exact",
        "Approximate equation": "approx",
        "Approximate epitope": "epitope",
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