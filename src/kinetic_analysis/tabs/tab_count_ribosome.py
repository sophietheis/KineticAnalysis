from dash import html, dcc
import dash_bootstrap_components as dbc

from kineticanalysis.tabs.tab_utils import title_h4, button, input_case_with_tooltip
from kineticanalysis.utils.texts import t


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([t("tab_ribosome.intro"), html.Br(),]),
            dcc.Markdown(t("tab_ribosome.math"), mathjax=True),
        ]),
        html.Br(),
        # Upload dataframes
        title_h4("Upload data"),
        html.Br(),

        dbc.Row([
            dbc.Col([
                html.H5("Single protein fluorescence", ),
                html.Div([
                    # Choose a csv file
                    html.Label(t("tab_ribosome.select_csv_single")),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_single_prot',
                                children=button(t("button.upload_csv")),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col([
                            dbc.Spinner(
                                children=[
                                    html.Div(id="loading_data_single_prot")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-single_prot'),
                ]),

                html.Br(),
                html.Label("Single protein fluorescence dataFrame visualisation"),
                html.Div(id="table_container_single", children=[]),
                html.Br(),
            ], width=6),

            dbc.Col([
                html.H5("Polysome fluorescence", ),
                html.Div([
                    # Choose a csv file
                    html.Label(t("tab_ribosome.select_csv_polysome")),
                    html.Br(),
                    dbc.Row([
                        dbc.Col([
                            dcc.Upload(
                                id='browse_directory_polysome',
                                children=button(t("button.upload_csv")),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col([
                            dbc.Spinner(
                                children=[
                                    html.Div(id="loading_data_polysome")],
                                size="lm",  # "sm"
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='selected-file-output-polysome'),
                ]),

                html.Br(),
                html.Label("Polysome fluorescence dataFrame visualisation"),
                html.Div(id="table_container_polysome", children=[]),

            ], width=6),
        ]),
        html.Br(),
        html.Br(),
        html.P(t("tab_ribosome.after_upload")),
        html.Br(),
        html.Br(),
        # Define parameter of the simulation
        dbc.Row([
            dbc.Col([input_case_with_tooltip(child="Protein length (aa)",
                                             input_id='param_prot_length_rib',
                                             tooltip_id="faq_param_prot_length_count_rib",
                                             tooltip_text=t("tooltips.prot_length"),
                                             input_type='number',
                                             input_value=t("default_input_value.prot_length")),
                    input_case_with_tooltip(child="Suntag length (aa)",
                                            input_id='param_suntag_length_rib',
                                            tooltip_id="faq_param_suntag_length",
                                            tooltip_text=t("tooltips.suntag_length"),
                                            input_type='number',
                                            input_value=t("default_input_value.suntag_length")),
                     ], width=3),


            # Calculate ribosome density
            dbc.Col(button(child=t("button.density"),
                           id_='btn_calculate_ribosome_density',
                           classname="mr-1"),
                    width=5)

        ]),
        html.Div(id='output_ribosome_density'),
        dcc.Download(id="download_csv"),
        dbc.Row([
            dbc.Col(children=[
                dcc.Graph(id='ribosome-plot'),
                ]),
        ]),
        dbc.Row([
            dbc.Col(children=[
                dcc.Graph(id='ribosome-density-plot'),
                ]),
            ]),
    ]),
    )
