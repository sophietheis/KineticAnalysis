from dash import html, dcc
import dash_bootstrap_components as dbc

from kineticanalysis.utils.texts import t
from kineticanalysis.tabs.tab_utils import button, input_case, title_h4


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([
                t("tab_msd.intro"),
                html.Br(),
                ]),
            dcc.Markdown(t("tab_msd.math"),
                         mathjax=True),
            html.P([
                t("tab_msd.msd_col_note"),
                html.Br(),
            ]),
        ]),
        html.Br(),
        # Upload dataframes
        title_h4("Upload data"),

        html.Br(),

        dbc.Row([
            dbc.Col([
                html.Div([
                    # Choose a csv file
                    html.Label("Choose your csv file. "),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_msd',
                                children=button(child=t("button.upload_csv")),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col([
                            dbc.Spinner(
                                children=[
                                    html.Div(id="loading_data_msd")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-msd'),
                ]),

                html.Br(),
                html.Div(id="table_container_msd", children=[]),
                html.Br(),
            ], width=6),

        ]),
        html.Br(),
        html.Br(),
        # CHOOSE COLUMN NAME
        title_h4("Select column name for the analysis"),
        html.Br(),
        dbc.Row([
            # Track ID column name
            dbc.Col(input_case(child="Track id",
                               input_id='col_ID',
                               input_type='text',
                               input_value=t("default_input_value.id"))),
            # X column name
            dbc.Col(input_case(child="x",
                               input_id='col_x',
                               input_type='text',
                               input_value=t("default_input_value.x"))),

            # Y column name
            dbc.Col(input_case(child="y",
                               input_id='col_y',
                               input_type='text',
                               input_value=t("default_input_value.y"))),

            # Z column name
            dbc.Col(input_case(child="z",
                               input_id='col_z',
                               input_type='text',
                               input_value=t("default_input_value.z"))),
        ]),

        html.Br(),
        html.Br(),
        dbc.Row([
            dbc.Col([
                button(child=t("button.msd"),
                       id_='btn_calculate_MSD',
                       classname="mr-2"),
            ], width=5)
        ]),
        html.Div(id='output_MSD'),
        dcc.Download(id="download_csv_msd"),
    ]),
    )
