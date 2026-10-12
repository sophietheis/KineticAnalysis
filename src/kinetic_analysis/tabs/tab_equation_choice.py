from dash import html, dcc
import dash_bootstrap_components as dbc

from kineticanalysis.utils.texts import t
from kineticanalysis.tabs.tab_utils import button, input_case_with_tooltip


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.P([t("tab_equation_choice.intro"),
                html.Br(),
                t("tab_equation_choice.reference"),
                html.Br(),
                ]),
        ]),
        dbc.Row([
            html.P([t("equations.intro_exact")]),
            dcc.Markdown(children=t("equations.g_exact"), mathjax=True),
        ]),
        dbc.Row([
            html.P([t("equations.intro_epitope")]),
            dcc.Markdown(children=t("equations.g_epitope"), mathjax=True),
            html.Br(),
            html.P([t("equations.intro_approx")]),
            dcc.Markdown(children=t("equations.g_approx"), mathjax=True),
        ]),
        html.Br(),
        html.Br(),
        # Define parameter of the simulation
        dbc.Row([
            dbc.Col([
                input_case_with_tooltip(child=t("labels.prot"),
                                        input_id='param_prot_length_choice',
                                        tooltip_id="faq_param_prot_length_choice",
                                        tooltip_text=t("tooltips.prot_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.prot_length")),

                input_case_with_tooltip(child=t("labels.suntag"),
                                        input_id='param_suntag_length_choice',
                                        tooltip_id="faq_param_suntag_length_choice",
                                        tooltip_text=t("tooltips.suntag_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.suntag_length")),

                input_case_with_tooltip(child=t("labels.nb_stemloop"),
                                        input_id='param_nb_suntag_choice',
                                        tooltip_id="faq_param_nb_suntag_choice",
                                        tooltip_text=t("tooltips.nb_stemloop"),
                                        input_type='number',
                                        input_value=t("default_input_value.nb_stemloop")),

                input_case_with_tooltip(child=t("labels.estimate_elongation_rate"),
                                        input_id='param_elongation_rate_choice',
                                        tooltip_id="faq_param_translation_rate_choice",
                                        tooltip_text=t("tooltips.elongation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.elongation_rate")),

                input_case_with_tooltip(child=t("labels.estimate_initiation_rate"),
                                        input_id='param_initiation_rate_choice',
                                        tooltip_id="faq_param_initiation_rate_choice",
                                        tooltip_text=t("tooltips.initiation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.initiation_rate")),

                input_case_with_tooltip(child=t("labels.tau"),
                                        input_id='param_tau_choice',
                                        tooltip_id="faq_param_dt_choice",
                                        tooltip_text=t("tooltips.tau"),
                                        input_type='number',
                                        input_value=t("default_input_value.tau")),

                html.Br(),
                html.Br(),
                dbc.Row([
                    button(t("button.show_contribution"), id_='show-contribution-btn')
                ]),
            ], width=4),

            # Show plot contribution
            dbc.Col([
                dcc.Graph(id='equation-plot'),

            ], width=8)
        ]),
        dbc.Row([
            dcc.Graph(id="equation-curve"),
        ]),
    ]),
    )
