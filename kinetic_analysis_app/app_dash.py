import os
import time
import threading
import subprocess

import numpy as np
import pandas as pd

from dash import Dash, html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_spinner

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from kinetic_function import (generate_profile,
                              generate_track,
                              single_track_analysis)

app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
app.title = "Kinetic analysis app"

# Global variables to store states
app.data = {
    'directory': None,
    'csv_files': [],
    'fig': None,
    'selected_file': None,
}


app.layout = dbc.Container([
    html.H1("Kinetic Analysis"),
    html.P("This tool is made to estimate the kinetic parameters of translation (initiation rate and elongation rate). "),
    html.P("You can find different tabs: "),
    html.Li("Track generator: generate tracks to test parameters"),
    html.Li("Track analysis: analyse tracks"),
    html.Br(),
    dbc.Tabs([
        dbc.Tab(label="Generate Tracks", children=[
            dbc.Row([
                html.Label("Choose Directory to Store Output"),
                html.Br(),
                dbc.Button("Select folder", id="add", className="mr-2", style={"width": "150px"},),
                # html.Div(id="output"),
                # dcc.Input(id='directory-input', type='text', placeholder='Directory path...'),
                # html.Button('Valid path', id='browse-dir-btn'),
                html.Div(id='directory-output', style={'margin-top': '10px'}),
                html.Br(),
                html.Br(),
                dbc.Col([
                    html.Div([
                        html.P("Protein length (aa)", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param1', type='number', value=490),
                    ]),
                    html.Div([
                        html.P("Suntag length (aa)", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param2', type='number', value=796),
                    ]),
                    html.Div([
                        html.P("Number of suntag", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param3', type='number', value=32),
                    ]),
                    html.Div([
                        html.P("Fluorescence one suntag", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param4', type='number', value=4),
                    ]),
                    html.Div([
                        html.P("Translation rate", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param5', type='number', value=24),
                    ]),
                    html.Div([
                        html.P("Binding rate", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param6', type='number', value=0.05),
                    ]),
                    html.Div([
                        html.P("Retention time", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param7', type='number', value=0),
                    ]),
                    html.Div([
                        html.P("Suntag position (begin or end)", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param8', type='text', value='begin'),
                    ]),
                    html.Div([
                        html.P("Number of tracks", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param9', type='number', value=100),
                    ]),
                    html.Div([
                        html.P("File name to save", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param10', type='text', value='datas'),
                    ]),
                    html.Br(),
                    # Generate Button and Spinner Side by Side
                    dbc.Row([
                        dbc.Col([
                            dbc.Button('Start Generate Tracks', id='start-gen-tracks-btn'),
                        ], width="auto"),

                        dbc.Col([
                            dbc.Spinner(
                                children=[html.Div(id="loading_generate")],
                                size="sm", color="primary", type="border", spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ], align="center", style={"margin-top": "10px"}),

                    html.Div(id='gen-tracks-output'),
                ], width=3),

                dbc.Col([
                    dbc.Button('Show Profile', id='show-profile-btn', className="mr-1"),
                    dcc.Graph(id='profile-plot'),

                ], width=5)
            ]),
        ]),


        dbc.Tab(label="Analyze Tracks", children=[
            html.Div([
                html.Label("Choose directory where your file is"),
                dcc.Input(id='directory-analyze-input', type='text', placeholder='Directory path...'),
                html.Button('Browse', id='browse-analyze-dir-btn'),
                html.Div(id='directory-analyze-output', style={'margin-top': '10px'}),
                html.Br(),
                dcc.Dropdown(id='file-dropdown', options=[], placeholder="Select a file..."),
                html.Button('Select File', id='select-file-btn'),
                html.Div(id='selected-file-output'),
                html.Div([
                    html.Label("dt"),
                    dcc.Input(id='dt-param', type='number', value=3),
                ]),
                html.Div([
                    html.Label("Protein length (aa)"),
                    dcc.Input(id='prot-length-param', type='number', value=800),
                ]),
                html.Div([
                    html.Label("File name to save"),
                    dcc.Input(id='save-results-name', type='text', value='datas_results'),
                ]),
                html.Br(),
                html.Button('Start Analyze Tracks', id='start-analyze-btn'),
                html.Div(id='analyze-output'),
            ]),
        ]),
    ]),
])

# Callbacks

# @app.callback(
#     Output('directory-output', 'children'),
#     Input('browse-dir-btn', 'n_clicks'),
#     State('directory-input', 'value')
# )
# def browse_directory(n_clicks, directory):
#     if n_clicks:
#         app.data['directory'] = directory
#         return f"Directory chosen: {directory}"
#     raise PreventUpdate

@app.callback(
    Output("directory-output", "children"),
    [Input("add", "n_clicks")],
)
def add(n_clicks):
    if n_clicks :
        command = ['python3', './local.py']
        p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out, err = p.communicate()
        app.data['directory'] = str(out.decode("utf-8"))[:-1] + "/"
        return f"Directory chosen: {app.data['directory']}"


@app.callback(
    Output('directory-analyze-output', 'children'),
    Input('browse-analyze-dir-btn', 'n_clicks'),
    State('directory-analyze-input', 'value')
)
def browse_directory_analyze(n_clicks, directory):
    if n_clicks:
        app.data['directory'] = directory
        load_csv_files(directory)
        return f"Directory chosen: {directory}"
    raise PreventUpdate


@app.callback(
    Output('file-dropdown', 'options'),
    Input('directory-analyze-output', 'children')
)
def load_csv_files(directory):
    if app.data['directory']:
        app.data['csv_files'] = [
            {'label': file, 'value': file}
            for file in os.listdir(app.data['directory']) if file.endswith('.csv')
        ]
        return app.data['csv_files']
    return []


@app.callback(
    Output('selected-file-output', 'children'),
    Input('select-file-btn', 'n_clicks'),
    State('file-dropdown', 'value')
)
def select_file(n_clicks, selected_file):
    if n_clicks and selected_file:
        app.data['selected_file'] = selected_file
        return f"You selected: {selected_file}"
    raise PreventUpdate


@app.callback(
    Output('profile-plot', 'figure'),
    Input('show-profile-btn', 'n_clicks'),
    State('param1', 'value'),
    State('param2', 'value'),
    State('param3', 'value'),
    State('param4', 'value'),
    State('param5', 'value'),
    State('param6', 'value'),
    State('param7', 'value'),
    State('param8', 'value'),
)
def update_profile_plot(n_clicks, *params):
    if n_clicks:
        try:
            # Generate profile
            x, y = generate_profile(*map(float, params[:7]), params[7])
            # Create the figure
            figure = make_subplots(rows=2, cols=1, subplot_titles=('One protein fluo profile',  'One track fluo profile'))
            figure.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Profile'),row=1, col=1)
            figure.update_xaxes(title_text='Time', row=1, col=1)
            figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

            if params[7] == "begin":
                suntag_pos = 0
            else:
                suntag_pos = -1
            x, y, _ = generate_track(float(params[0]),
                                     float(params[1]),
                                     float(params[2]),
                                     float(params[3]),
                                     float(params[4]),
                                     float(params[5]),
                                     float(params[6]),
                                     suntag_pos)

            figure.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Profile'), row=2, col=1)
            figure.update_xaxes(title_text='Time', row=2, col=1)
            figure.update_yaxes(title_text='Fluorescence', row=2, col=1)

            figure.update_layout(width=1000, height=800,)
            return figure
        except Exception as e:
            return {
                'data': [],
                'layout': go.Layout(title='Error', xaxis={'title': 'Time'}, yaxis={'title': 'Fluorescence'})
            }
    raise PreventUpdate


@app.callback(
    Output('gen-tracks-output', 'children'),
    Output('loading_generate', 'children'),
    Input('start-gen-tracks-btn', 'n_clicks'),
    State('param1', 'value'),
    State('param2', 'value'),
    State('param3', 'value'),
    State('param4', 'value'),
    State('param5', 'value'),
    State('param6', 'value'),
    State('param7', 'value'),
    State('param8', 'value'),
    State('param9', 'value'),
    State('param10', 'value'),
)
def start_generate_tracks(n_clicks, *params):
    if n_clicks:
        try:
            # Generate all tracks and save it
            first_time = True
            if params[7] == "begin":
                suntag_pos = 0
            else:
                suntag_pos = -1

            for i in range(int(params[8])):
                x_global, y_global, y_start_prot = generate_track(float(params[0]),
                                                                 float(params[1]),
                                                                 float(params[2]),
                                                                 float(params[3]),
                                                                 float(params[4]),
                                                                 float(params[5]),
                                                                 float(params[6]),
                                                                 suntag_pos)
                if first_time:
                    datas = pd.DataFrame({"FRAME": x_global,
                                          "MEAN_INTENSITY_CH1": y_global,
                                          "TRACK_ID": i,
                                          "RETENTION_TIME": float(params[6]),
                                          })
                    first_time = False
                else:
                    datas = pd.concat([datas,
                                       pd.DataFrame({"FRAME": x_global,
                                                     "MEAN_INTENSITY_CH1": y_global,
                                                     "TRACK_ID": i,
                                                     "RETENTION_TIME": float(params[6]),
                                                     })], ignore_index=True)

            print(app.data["directory"])
            datas.to_csv(os.path.join(app.data['directory'],  params[9] + ".csv"))

            return  "Tracks generated and saved successfully!", None
        except Exception as e:
            return  f"Error: {str(e)}", None
    raise PreventUpdate


@app.callback(
    Output('analyze-output', 'children'),
    Input('start-analyze-btn', 'n_clicks'),
    State('dt-param', 'value'),
    State('prot-length-param', 'value'),
)
def start_analyze_tracks(n_clicks, dt, prot_length):
    if n_clicks:
        try:
            # Implement track analysis here
            # Simulate analysis for now
            time.sleep(2)  # Simulate time taken for processing
            return "Analysis completed and saved successfully!"
        except Exception as e:
            return f"Error: {str(e)}"
    raise PreventUpdate


if __name__ == '__main__':
    app.run_server(debug=True, port=8050)

