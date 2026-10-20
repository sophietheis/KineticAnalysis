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
                                id='countrib_browse_single_prot',
                                children=button(t("button.upload_csv")),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col([
                            dbc.Spinner(
                                children=[
                                    html.Div(id="countrib_load_single_prot")],
                                size="sm",
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='countrib_select_single_prot'),
                ]),

                html.Br(),
                html.Label("Single protein fluorescence dataFrame visualisation"),
                html.Div(id="countrib_table_container_single_prot", children=[]),
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
                                id='countrib_browse_polysome',
                                children=button(t("button.upload_csv")),
                                multiple=False,
                            ),
                        ], width="auto"),
                        dbc.Col([
                            dbc.Spinner(
                                children=[
                                    html.Div(id="countrib_load_polysome")],
                                size="lm",  # "sm"
                                color="primary",
                                type="border",
                                spinner_style={"margin-left": "10px"}
                            )
                        ], width="auto"),
                    ]),
                    html.Div(id='countrib_select_polysome'),
                ]),

                html.Br(),
                html.Label("Polysome fluorescence dataFrame visualisation"),
                html.Div(id="countrib_table_container_polysome", children=[]),

            ], width=6),
        ]),
        html.Br(),
        html.Br(),
        html.P(t("tab_ribosome.after_upload")),
        html.Br(),
        html.Br(),
        # Define parameter of the simulation
        dbc.Row([
            dbc.Col([input_case_with_tooltip(child=t("labels.prot"),
                                             input_id='countrib_prot_length',
                                             tooltip_id="countrib_faq_prot_length",
                                             tooltip_text=t("tooltips.prot_length"),
                                             input_type='number',
                                             input_value=t("default_input_value.prot_length")),
                    input_case_with_tooltip(child=t("labels.epitope"),
                                            input_id='countrib_epitope_length',
                                            tooltip_id="countrib_faq_epitope_length",
                                            tooltip_text=t("tooltips.suntag_length"),
                                            input_type='number',
                                            input_value=t("default_input_value.suntag_length")),
                     ], width=3),


            # Calculate ribosome density
            dbc.Col(button(child=t("button.density"),
                           id_='countrib_button_calculate',
                           classname="mr-1"),
                    width=5)

        ]),
        html.Div(id='countrib_output'),
        dcc.Download(id="countrib_download_csv"),
        dbc.Row([
            dbc.Col(children=[
                dcc.Graph(id='countrib_display_plot'),
                ]),
        ]),
        dbc.Row([
            dbc.Col(children=[
                dcc.Graph(id='countrib_display_density_plot'),
                ]),
            ]),
    ]),
    )
