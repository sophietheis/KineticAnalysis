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
                                        input_id='eqchoice_prot_length',
                                        tooltip_id="eqchoice_faq_prot_length",
                                        tooltip_text=t("tooltips.prot_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.prot_length")),

                input_case_with_tooltip(child=t("labels.epitope"),
                                        input_id='eqchoice_epitope_length',
                                        tooltip_id="eqchoice_faq_epitope_length",
                                        tooltip_text=t("tooltips.suntag_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.suntag_length")),

                input_case_with_tooltip(child=t("labels.nb_stemloop"),
                                        input_id='eqchoice_epitope_repeat',
                                        tooltip_id="eqchoice_faq_epitope_repeat",
                                        tooltip_text=t("tooltips.nb_stemloop"),
                                        input_type='number',
                                        input_value=t("default_input_value.nb_stemloop")),

                input_case_with_tooltip(child=t("labels.estimate_elongation_rate"),
                                        input_id='eqchoice_elongation_rate',
                                        tooltip_id="eqchoice_faq_elongation_rate",
                                        tooltip_text=t("tooltips.elongation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.elongation_rate")),

                input_case_with_tooltip(child=t("labels.estimate_initiation_rate"),
                                        input_id='eqchoice_initiation_rate',
                                        tooltip_id="eqchoice_faq_initiation_rate",
                                        tooltip_text=t("tooltips.initiation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.initiation_rate")),

                input_case_with_tooltip(child=t("labels.tau"),
                                        input_id='eqchoice_tau',
                                        tooltip_id="eqchoice_faq_tau",
                                        tooltip_text=t("tooltips.tau"),
                                        input_type='number',
                                        input_value=t("default_input_value.tau")),

                html.Br(),
                html.Br(),
                dbc.Row([
                    button(t("button.show_contribution"), id_='eqchoice_button_display')
                ]),
            ], width=4),

            # Show plot contribution
            dbc.Col([
                dcc.Graph(id='eqchoice_plot'),

            ], width=8)
        ]),
        dbc.Row([
            dcc.Graph(id="eqchoice_curve"),
        ]),
    ]),
    )
