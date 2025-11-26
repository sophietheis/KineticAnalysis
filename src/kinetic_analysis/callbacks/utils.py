from dash import dash_table

def generate_table(dataframe, max_rows=10):
    table = dash_table.DataTable(
        data=dataframe.to_dict('records'),
        columns=[{"name": i, "id": i} for i in
                 dataframe.columns],
        page_size=max_rows,
        style_table={'width': '800px', 'overflowX': 'auto'},
    )
    return table
