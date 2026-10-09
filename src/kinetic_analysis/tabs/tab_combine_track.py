from dash import html, dcc
import dash_bootstrap_components as dbc

from kineticanalysis.utils.texts import t
from kineticanalysis.tabs.tab_utils import button, input_case, title_h4


def layout():
    return (
        html.Br(),
        # FILE IMPORT AND DISPLAY
        title_h4("Upload data"),

        html.Br(),
        dbc.Row([
            dbc.Col(children=[
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file to analyse "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col(children=[
                            dcc.Upload(
                                id='browse_directory_combine',
                                children=button(t("button.upload_csv")),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col(children=[
                            dbc.Spinner(
                                children=[
                                    html.Div(id="loading_data_combine")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-combine'),
                ]),
            ], width=3),

            dbc.Col(children=[
                html.Label("DataFrame Visualisation"),
                html.Div(id="table-container_combine", children=[]),
            ], width=9)
        ]),

        html.Br(),

        dbc.Row([
            html.P(children=["There is : ",
                             html.Span(id='nb_tracks_init', children=''),
                             " track(s)."],
                   className="mb-0"),
        ]),

        html.Br(),
        title_h4("Parameters for combining tracks"),

        dbc.Row([
            # nb track in the combined track
            dbc.Col(input_case(child="How many tracks in the combined tracks?",
                               input_id='nb_tracks',
                               input_type='number',
                               input_value="2")),

            # number of new tracks created
            dbc.Col(input_case(child="How many new track created ?",
                               input_id='nb_new_tracks',
                               input_type='number',
                               input_value="10")),
        ]),

        html.Br(),

        dbc.Row([
            html.Div([
                dbc.Col(children=[
                    dcc.Store(id="start2", data=""),
                    dcc.Store(id="complete2", data=""),
                    dbc.Button(dbc.Spinner(
                        html.Span(children="Combine tracks",
                                  id="loading_combination")),
                        id="combine-tracks-btn"),
                ], width=3),
            ]),
            html.Div(id='combine-tracks-output'),
            dcc.Download(id="download-csv3"),
        ]),

    )
