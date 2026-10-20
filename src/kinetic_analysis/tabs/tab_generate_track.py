from dash import html, dcc
import dash_bootstrap_components as dbc

from kineticanalysis.tabs.tab_utils import button, input_case_with_tooltip, input_case
from kineticanalysis.utils.texts import t


def layout():
    return (html.Div([
        # Explanation at the beginning of the page
        dbc.Row([
            html.Br(),
            html.P([t("tab_generate_track.intro"),
                    html.Br(),
                    t("tab_generate_track.intro2"),
                    ]),
            html.Br(),
        ]),
        dbc.Row([
            html.P([t("tab_generate_track.eq_fluo")]),
            dcc.Markdown(children=t("tab_generate_track.eq_fluo_markdown"), mathjax=True),
        ]),
        dbc.Row([
            html.P([t("tab_generate_track.eq_total_fluo")]),
            dcc.Markdown(children=t("tab_generate_track.eq_total_fluo_markdown"), mathjax=True),

        ]),

        # Define parameter of the simulation
        dbc.Row([
            dbc.Col(children=[
                input_case_with_tooltip(child=t("labels.prot"),
                                        input_id='gentrack_prot_length',
                                        tooltip_id="gentrack_faq_prot_length",
                                        tooltip_text=t("tooltips.prot_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.prot_length")),

                input_case_with_tooltip(child=t("labels.epitope"),
                                        input_id='gentrack_epitope_length',
                                        tooltip_id="gentrack_faq_epitope_length",
                                        tooltip_text=t("tooltips.suntag_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.suntag_length")),

                input_case_with_tooltip(child=t("labels.nb_stemloop"),
                                        input_id='gentrack_epitope_repeat',
                                        tooltip_id="gentrack_faq_epitope_repeat",
                                        tooltip_text=t("tooltips.nb_stemloop"),
                                        input_type='number',
                                        input_value=t("default_input_value.nb_stemloop")),

                input_case_with_tooltip(child=t("labels.fluo_one_suntag"),
                                        input_id='gentrack_epitope_fluo',
                                        tooltip_id="gentrack_faq_epitope_fluo",
                                        tooltip_text=t("tooltips.nb_stemloop"),
                                        input_type='number',
                                        input_value=t("default_input_value.fluo_one_suntag")),

                input_case_with_tooltip(child=t("labels.elongation_rate"),
                                        input_id='gentrack_elongation_rate',
                                        tooltip_id="gentrack_faq_elongation_rate",
                                        tooltip_text=t("tooltips.elongation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.elongation_rate")),

                input_case_with_tooltip(child=t("labels.initiation_rate"),
                                        input_id='gentrack_initiation_rate',
                                        tooltip_id="gentrack_faq_initiation_rate",
                                        tooltip_text=t("tooltips.initiation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.initiation_rate")),

                input_case_with_tooltip(child=t("labels.footprint"),
                                        input_id='gentrack_footprint',
                                        tooltip_id="gentrack_faq_footprint",
                                        tooltip_text=t("tooltips.footprint"),
                                        input_type='number',
                                        input_value=t("default_input_value.footprint")),

                input_case_with_tooltip(child=t("labels.retention_time"),
                                        input_id='gentrack_retention_time',
                                        tooltip_id="gentrack_faq_retention_time",
                                        tooltip_text=t("tooltips.retention_time"),
                                        input_type='number',
                                        input_value=t("default_input_value.retention_time")),

                html.Div([
                    html.P(children=[t("labels.pos_suntag"),
                                     html.Span(className="fas fa-question-circle",
                                               id="gentrack_faq_epitope_position",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children=t("tooltips.pos_suntag"),
                                target="gentrack_faq_epitope_position"),
                    dcc.Dropdown(id='gentrack_epitope_position',
                                 options=[
                                     {'label': 'Begin', 'value': 'begin'},
                                     {'label': 'End', 'value': 'end'}],
                                 value='begin',  # Default value
                                 clearable=False,
                                 # Prevents the user from clearing the selection
                                 style={'width': '200px'}
                                 # Adjust width if needed
                                 ),
                ]),
                input_case_with_tooltip(child=t("labels.noise"),
                                        input_id='gentrack_noise',
                                        tooltip_id="gentrack_faq_noise",
                                        tooltip_text=t("tooltips.noise"),
                                        input_type='number',
                                        input_value=t("default_input_value.noise")),

                input_case_with_tooltip(child=t("labels.dt"),
                                        input_id='gentrack_dt',
                                        tooltip_id="gentrack_faq_dt",
                                        tooltip_text=t("tooltips.dt"),
                                        input_type='number',
                                        input_value=t("default_input_value.dt")),

                input_case_with_tooltip(child=t("labels.track_length"),
                                        input_id='gentrack_track_length',
                                        tooltip_id="gentrack_faq_track_length",
                                        tooltip_text=t("tooltips.track_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.track_length")),

                input_case_with_tooltip(child=t("labels.nb_tracks"),
                                        input_id='gentrack_nb_track',
                                        tooltip_id="gentrack_faq_nb_track",
                                        tooltip_text=t("tooltips.nb_tracks"),
                                        input_type='number',
                                        input_value=t("default_input_value.nb_tracks")),

                input_case(child=t("labels.filename"),
                           input_id='gentrack_filename',
                           input_type='text',
                           input_value=t("default_input_value.filename_data")),

            ], width=4),

            # Graph area
            dbc.Col(children=[
                dcc.Graph(id='gentrack_plot_profile'),
            ], width=8),
        ]),

        html.Br(),
        html.Br(),

        # Button row
        dbc.Row([
            # Show plot profile button
            dbc.Col(children=[button(child=t("button.show_profile"),
                                     id_='gentrack_button_profile')],
                    width=2),

            # Generate track Button and Spinner Side by Side
            dbc.Col([
                html.Div([
                    dbc.Col(children=[
                        dcc.Store(id="gentrack_start", data=""),
                        dcc.Store(id="gentrack_complete", data=""),
                        dbc.Button(dbc.Spinner(
                            html.Span(children="Generate tracks",
                                      id="gentrack_loading")),
                            id="gentrack_button_start"),
                    ], width=3),
                ]),
                html.Div(id='gentrack_output'),
                dcc.Download(id="gentrack_download_csv"),
            ]),
        ]),
    ]),
    )
