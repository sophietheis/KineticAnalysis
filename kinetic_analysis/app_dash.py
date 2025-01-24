import os
import time
import threading
import webview

from threading import Thread

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog

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
    'directory_generation': None,
    'directory_analysis': None,
    'csv_files': [],
    'fig': None,
    'selected_file': None,
}


app.layout = dbc.Container([
    html.H1("Kinetic Analysis"),
    html.P("This tool is made to estimate the kinetic parameters of translation (initiation rate and elongation rate). "),
    html.P("You can find different tabs: "),
    html.Li("Track generator: generate tracks to test parameters"),
    html.Li("Track analysis simulation : analyse tracks simulation"),
    html.Li("Track analysis in vivo : analyse tracks from in vivo "
            "experiments"),
    html.Br(),
    dbc.Tabs([

        # Generate track tab
        dbc.Tab(label="Generate Tracks", children=[
            dbc.Row([
                html.Label("Choose Directory to Store Output"),
                html.Br(),
                dbc.Button("Select folder", id="add", className="mr-2", style={"width": "150px"},),
                html.Div(id='directory-output', style={'margin-top': '10px'}),
                html.Br(),
                html.Br(),
                dbc.Col([
                    html.Div([
                        html.P("Protein length (aa)", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param1', type='number', value=490, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("Suntag length (aa)", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param2', type='number', value=796, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("Number of suntag", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param3', type='number', value=32, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("Fluorescence one suntag", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param4', type='number', value=4, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("Translation rate (aa/sec)", style={"height":
                                                                "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param5', type='number', value=24, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("Initiation rate (ribosome/sec)", style={
                            "height":
                                                                   "auto",
                                                       "margin-bottom": "auto"}),
                        dcc.Input(id='param6', type='number', value=1, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("Retention time (sec)", style={"height":
                                                                  "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param7', type='number', value=0, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("Suntag position (begin or end)", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Dropdown(id='param8',
                                     options=[{'label': 'Begin', 'value': 'begin'},
                                              {'label': 'End', 'value': 'end'}],
                                     value='begin',  # Default value
                                     clearable=False,  # Prevents the user from clearing the selection
                                     style={'width': '200px'}  # Adjust width if needed
                        ),
                    ]),
                    html.Div([
                        html.P("Number of tracks", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param9', type='number', value=100, style={'width': '200px'}),
                    ]),
                    html.Div([
                        html.P("File name to save", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='param10', type='text', value='datas', style={'width': '200px'}),
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

        # Analyse track simulation tab
        dbc.Tab(label="Analyze Tracks Simulation", children=[
            html.Div([

                # Choose a directory
                html.Label("Choose Directory where your file is"),
                html.Br(),
                dbc.Button("Select folder", id="browse_directory_analyze", className="mr-2", style={"width": "150px"},),
                html.Div(id='directory-analyze-output', style={'margin-top': '10px'}),
                html.Br(),
                # Select a file inside the directory
                dcc.Dropdown(id='file-dropdown', options=[], placeholder="Select a file...", style={"width": "150px"},),
                dbc.Button('Validate file', id='select-file-btn', className="mr-2", style={"width": "150px"},),
                html.Div(id='selected-file-output'),
                html.Br(),
                html.Br(),
                html.Br(),
                dbc.Col([
                    html.Div([
                        html.P("dt", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='dt-param', type='number', value=3),
                    ]),
                    html.Div([
                        html.P("Protein length (aa)", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='prot-length-param', type='number', value=800),
                    ]),
                    html.Div([
                        html.P("File name to save", style={"height": "auto", "margin-bottom": "auto"}),
                        dcc.Input(id='save-results-name', type='text', value='datas_results'),
                    ]),
                    html.Br(),
                    # Generate Button and Spinner Side by Side
                    dbc.Row([
                        dbc.Col([
                            dbc.Button('Start Analyze Tracks', id='start-analyze-btn', className="mr-2", style={"width": "150px"},),
                        ], width="auto"),

                        dbc.Col([
                            dbc.Spinner(
                                children=[html.Div(id="loading_analysis")],
                                size="sm", color="primary", type="border", spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                        html.Div(id='analyze-output'),
                    ], align="center", style={"margin-top": "10px"}),
                ], width=3),
            ]),
        ]),

        # Analyse track in vivo tab
        dbc.Tab(label="Analyze Tracks in vivo", children=[
            html.Div([

                # Choose a directory
                html.Label("Choose Directory where your file is"),
                html.Br(),
                dbc.Button("Select folder",
                           id="browse_directory_analyze_vivo",
                           className="mr-2", style={"width": "150px"},),
                html.Div(id='directory-analyze-output-vivo', style={
                    'margin-top': '10px'}),
                html.Br(),
                # Select a file inside the directory
                dcc.Dropdown(id='file-dropdown-vivo', options=[],
                             placeholder="Select a file...", style={"width": "150px"},),
                dbc.Button('Validate file', id='select-file-btn-vivo',
                           className="mr-2", style={"width": "150px"},),
                html.Div(id='selected-file-output-vivo'),
                html.Br(),
                html.Br(),
                html.Br(),
        #         dbc.Col([
        #             html.Div([
        #                 html.P("dt", style={"height": "auto", "margin-bottom": "auto"}),
        #                 dcc.Input(id='dt-param', type='number', value=3),
        #             ]),
        #             html.Div([
        #                 html.P("Protein length (aa)", style={"height": "auto", "margin-bottom": "auto"}),
        #                 dcc.Input(id='prot-length-param', type='number', value=800),
        #             ]),
        #             html.Div([
        #                 html.P("File name to save", style={"height": "auto", "margin-bottom": "auto"}),
        #                 dcc.Input(id='save-results-name', type='text', value='datas_results'),
        #             ]),
        #             html.Br(),
        #             # Generate Button and Spinner Side by Side
        #             dbc.Row([
        #                 dbc.Col([
        #                     dbc.Button('Start Analyze Tracks', id='start-analyze-btn', className="mr-2", style={"width": "150px"},),
        #                 ], width="auto"),
        #
        #                 dbc.Col([
        #                     dbc.Spinner(
        #                         children=[html.Div(id="loading_analysis")],
        #                         size="sm", color="primary", type="border", spinner_style={"margin-left": "10px"}
        #                     )
        #                 ], width="auto"),
        #                 html.Div(id='analyze-output'),
        #             ], align="center", style={"margin-top": "10px"}),
        #         ], width=3),
            ]),
        ]),
    ]),
])

# Callbacks


@app.callback(
    Output("directory-output", "children"),
    [Input("add", "n_clicks")],
)
def add(n_clicks):
    if n_clicks :
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder_selected = filedialog.askdirectory()
        root.destroy()
        if folder_selected:
            print(folder_selected)
        else:
            print(None)
        app.data['directory_generation'] = folder_selected
        return f"Directory chosen: {app.data['directory_generation']}"

@app.callback(
    Output('directory-analyze-output', 'children'),
    Input('browse_directory_analyze', 'n_clicks'),
)
def browse_directory_analyze(n_clicks):
    if n_clicks:
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder_selected = filedialog.askdirectory()
        root.destroy()
        if folder_selected:
            print(folder_selected)
        else:
            print(None)
        app.data['directory_analysis'] = folder_selected
        return f"Directory chosen: {app.data['directory_analysis']}"

@app.callback(
    Output('file-dropdown', 'options'),
    Input('directory-analyze-output', 'children')
)
def load_csv_files(directory):
    if app.data['directory_analysis']:
        app.data['csv_files'] = [
            {'label': file, 'value': file}
            for file in os.listdir(app.data['directory_analysis']) if file.endswith('.csv')
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
    Output('directory-analyze-output-vivo', 'children'),
    Input('browse_directory_analyze_vivo', 'n_clicks'),
)
def browse_directory_analyze_vivo(n_clicks):
    if n_clicks:
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder_selected = filedialog.askdirectory()
        root.destroy()
        if folder_selected:
            print(folder_selected)
        else:
            print(None)
        app.data['directory_analysis_vivo'] = folder_selected
        return f"Directory chosen: {app.data['directory_analysis_vivo']}"

@app.callback(
    Output('file-dropdown-vivo', 'options'),
    Input('directory-analyze-output-vivo', 'children')
)
def load_csv_files(directory):
    if app.data['directory_analysis_vivo']:
        app.data['csv_files'] = [
            {'label': file, 'value': file}
            for file in os.listdir(app.data['directory_analysis_vivo']) if file.endswith('.csv')
        ]
        return app.data['csv_files']
    return []

@app.callback(
    Output('selected-file-output-vivo', 'children'),
    Input('select-file-btn-vivo', 'n_clicks'),
    State('file-dropdown-vivo', 'value')
)
def select_file(n_clicks, selected_file):
    if n_clicks and selected_file:
        app.data['selected_file_vivo'] = selected_file
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
            x,y = generate_profile(float(params[0]),
                                     float(params[1]),
                                     float(params[2]),
                                     float(params[3]),
                                     float(params[4]),
                                     float(params[6]),
                                     params[7])
            # Create the figure
            figure = make_subplots(rows=3,
                                   cols=1,
                                   subplot_titles=('One protein fluo profile',
                                                   'One track fluo profile',
                                                   'Number of translation'))
            figure.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Profile'),row=1, col=1)
            figure.update_xaxes(title_text='Time', row=1, col=1)
            figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

            x, y, y_number = generate_track(float(params[0]),
                                     float(params[1]),
                                     float(params[2]),
                                     float(params[3]),
                                     float(params[4]),
                                     float(params[5]),
                                     float(params[6]),
                                     params[7])

            figure.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Profile'), row=2, col=1)
            figure.update_xaxes(title_text='Time', row=2, col=1)
            figure.update_yaxes(title_text='Fluorescence', row=2, col=1)

            figure.add_trace(go.Scatter(x=x, y=y_number, mode='lines',
                                        name='Profile'), row=3, col=1)
            figure.update_xaxes(title_text='Time', row=3, col=1)
            figure.update_yaxes(title_text='Number of translation', row=3,
                                col=1)
            figure.update_layout(width=1000, height=800,)
            return figure
        except Exception as e:
            print(e)
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

            for i in range(int(params[8])):
                x_global, y_global, y_start_prot = generate_track(float(params[0]),
                                                                 float(params[1]),
                                                                 float(params[2]),
                                                                 float(params[3]),
                                                                 float(params[4]),
                                                                 float(params[5]),
                                                                 float(params[6]),
                                                                 params[7])
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

            datas.to_csv(os.path.join(app.data['directory_generation'],  params[9] + ".csv"))

            return  "Tracks generated and saved successfully!", None
        except Exception as e:
            return  f"Error: {str(e)}", None
    raise PreventUpdate


@app.callback(
    Output('analyze-output', 'children'),
    Output('loading_analysis', 'children'),
    Input('start-analyze-btn', 'n_clicks'),
    State('file-dropdown', 'value'),
    State('dt-param', 'value'),
    State('prot-length-param', 'value'),
    State('save-results-name', 'value'),
)
def start_analyze_tracks(n_clicks, filename, *params):
    if n_clicks:
        try:
            # Read csv file
            datas = pd.read_csv(os.path.join(app.data['directory_analysis'], filename),
                                index_col="Unnamed: 0")
            dt = float(params[0])
            t = dt / 0.1
            prot_length = float(params[1])
            nb_track = len(np.unique(datas["TRACK_ID"]))

            first_time = True
            # Analyse all tracks and save it
            for i in range(nb_track):
                datas2 = datas[(datas["TRACK_ID"] == i)][::int(t)]

                (x,
                 y,
                 x_auto,
                 y_auto,
                 elongation_r,
                 translation_init_r,
                 perr) = single_track_analysis(datas2,
                                               i,
                                               delta_t=dt,
                                               protein_size=prot_length,
                                               normalise_intensity=1,
                                               normalize_auto=True,
                                               mm=None,
                                               lowpass_=False,
                                               cutoff=100,
                                               rtol=1e-1,
                                               method="linear",
                                               force_analysis=True,
                                               first_dot=True,
                                               simulation=True)
                if first_time:
                    results = pd.DataFrame({"elongation_r": elongation_r,
                                            "init_translation_r": translation_init_r,
                                            "dt": dt,
                                            "id": i, },
                                           index=[0])
                    first_time = False

                else:
                    results = pd.concat([results,
                                         pd.DataFrame({"elongation_r": elongation_r,
                                                       "init_translation_r": translation_init_r,
                                                       "dt": dt,
                                                       "id": i, }, index=[0])
                                         ], ignore_index=True)


            results.to_csv(os.path.join(app.data['directory_analysis'], params[2] + ".csv"))

            return "Analysis completed and saved successfully!", None
        except Exception as e:
            return f"Error: {str(e)}", None
    raise PreventUpdate


if __name__ == '__main__':
    # app.run_server(debug=True, port=8080)
    def run_app():
        app.run_server(debug=False, port=8080)
    t = Thread(target=run_app)
    t.daemon = True
    t.start()

    window = webview.create_window('Kinetic analysis', 'http://127.0.0.1:8080/',
                                   width=1800, height=1000,
)
    webview.start(debug=False)

