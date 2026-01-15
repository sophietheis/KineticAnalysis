from dash import dash_table


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
