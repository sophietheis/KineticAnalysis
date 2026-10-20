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
                                id='combine_browse_directoy',
                                children=button(t("button.upload_csv")),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col(children=[
                            dbc.Spinner(
                                children=[
                                    html.Div(id="combine_loading_data")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='combine_output_file'),
                ]),
            ], width=3),

            dbc.Col(children=[
                html.Label("DataFrame Visualisation"),
                html.Div(id="combine_table_container", children=[]),
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
                               input_id='combine_nb_tracks',
                               input_type='number',
                               input_value="2")),

            # number of new tracks created
            dbc.Col(input_case(child="How many new track created ?",
                               input_id='combine_nb_new_tracks',
                               input_type='number',
                               input_value="10")),
        ]),

        html.Br(),

        dbc.Row([
            html.Div([
                dbc.Col(children=[
                    dcc.Store(id="combine_start", data=""),
                    dcc.Store(id="combine_complete", data=""),
                    dbc.Button(dbc.Spinner(
                        html.Span(children="Combine tracks",
                                  id="combine_loading")),
                        id="combine_button_combine"),
                ], width=3),
            ]),
            html.Div(id='combine_output_tracks'),
            dcc.Download(id="combine_download_csv"),
        ]),

    )
