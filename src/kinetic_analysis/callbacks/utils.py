from dash import dash_table

import plotly.graph_objs as go


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

def resolve_solver_method(app):
    return {
        "Exact equation": "exact",
        "Approximate equation": "approx",
        "Approximate epitope": "epitope",
    }.get(app.data.get("solver"), "exact")


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