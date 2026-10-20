from dash import html
import dash_bootstrap_components as dbc

from kineticanalysis.tabs.tab_utils import title_h4, input_case, button
from kineticanalysis.utils.texts import t


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([t("tab_acquisition.intro"),
                    html.Br(),
                    t("tab_acquisition.intro1"),
                    html.Br(),
                    t("tab_acquisition.intro2"),
                    ]),
        ]),
        html.Br(),

        dbc.Row([
            dbc.Col([
                title_h4("Input value"),
            ], width=5),
            dbc.Col([
                title_h4("Output"),
            ], width=5),
        ]),
        html.Br(),

        dbc.Row([
            dbc.Col([
                input_case(child=t("labels.prot"),
                           input_id="param_prot_length_acquisition",
                           input_type='number',
                           input_value=t("default_input_value.prot_length")),

                input_case(child=t("labels.epitope"),
                           input_id="param_suntag_length_acquisition",
                           input_type='number',
                           input_value=t("default_input_value.suntag_length")),

                input_case(child=t("labels.nb_stemloop"),
                           input_id="param_stem_loops_acquisition",
                           input_type='number',
                           input_value=t("default_input_value.nb_stemloop")),

                html.Br(),

                input_case(child=t("labels.elongation_rate"),
                           input_id="param_translation_rate_acquisition",
                           input_type='number',
                           input_value=t("default_input_value.elongation_rate")),

                input_case(child=t("labels.nb_prot"),
                           input_id="param_full_translated_protein_acquisition",
                           input_type='number',
                           input_value=t("default_input_value.nb_prot")),

                html.Br(),
                button(child=t("button.acquisition"),
                       id_='btn_calculate_acquisition'),
            ], width=5),

            dbc.Col([
                html.Div(id='output_acquisition'),
            ])
        ]),

    ]),
    )
