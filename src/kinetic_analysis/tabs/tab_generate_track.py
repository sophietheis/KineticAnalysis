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
                                        input_id='param_prot_length_gen_track',
                                        tooltip_id="faq_param_prot_length_gen_track",
                                        tooltip_text=t("tooltips.prot_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.prot_length")),

                input_case_with_tooltip(child=t("labels.epitope"),
                                        input_id='param_suntag_length_gen_track',
                                        tooltip_id="faq_param_suntag_length_gen_track",
                                        tooltip_text=t("tooltips.suntag_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.suntag_length")),

                input_case_with_tooltip(child=t("labels.nb_stemloop"),
                                        input_id='param_nb_suntag_gen_track',
                                        tooltip_id="faq_param_nb_suntag_gen_track",
                                        tooltip_text=t("tooltips.nb_stemloop"),
                                        input_type='number',
                                        input_value=t("default_input_value.nb_stemloop")),

                input_case_with_tooltip(child=t("labels.fluo_one_suntag"),
                                        input_id='param_fluo_one_suntag_gen_track',
                                        tooltip_id="faq_param_fluo_one_suntag_gen_track",
                                        tooltip_text=t("tooltips.nb_stemloop"),
                                        input_type='number',
                                        input_value=t("default_input_value.fluo_one_suntag")),

                input_case_with_tooltip(child=t("labels.elongation_rate"),
                                        input_id='param_translation_rate_gen_track',
                                        tooltip_id="faq_param_translation_rate_gen_track",
                                        tooltip_text=t("tooltips.elongation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.elongation_rate")),

                input_case_with_tooltip(child=t("labels.initiation_rate"),
                                        input_id='param_initiation_rate_gen_track',
                                        tooltip_id="faq_param_initiation_rate_gen_track",
                                        tooltip_text=t("tooltips.initiation_rate"),
                                        input_type='number',
                                        input_value=t("default_input_value.initiation_rate")),

                input_case_with_tooltip(child=t("labels.footprint"),
                                        input_id='param_footprint_gen_track',
                                        tooltip_id="faq_param_footprint_gen_track",
                                        tooltip_text=t("tooltips.footprint"),
                                        input_type='number',
                                        input_value=t("default_input_value.footprint")),

                input_case_with_tooltip(child=t("labels.retention_time"),
                                        input_id='param_retention_time_gen_track',
                                        tooltip_id="faq_param_retention_time_gen_track",
                                        tooltip_text=t("tooltips.retention_time"),
                                        input_type='number',
                                        input_value=t("default_input_value.retention_time")),

                html.Div([
                    html.P(children=[t("labels.pos_suntag"),
                                     html.Span(className="fas fa-question-circle",
                                               id="faq_param_pos_suntag_gen_track",
                                               style={"cursor": "pointer",
                                                      "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(children=t("tooltips.pos_suntag"),
                                target="faq_param_pos_suntag_gen_track"),
                    dcc.Dropdown(id='param_pos_suntag_gen_track',
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
                                        input_id='param_noise_gen_track',
                                        tooltip_id="faq_param_noise_gen_track",
                                        tooltip_text=t("tooltips.noise"),
                                        input_type='number',
                                        input_value=t("default_input_value.noise")),

                input_case_with_tooltip(child=t("labels.dt"),
                                        input_id='param_dt_gen_track',
                                        tooltip_id="faq_param_dt_gen_track",
                                        tooltip_text=t("tooltips.dt"),
                                        input_type='number',
                                        input_value=t("default_input_value.dt")),

                input_case_with_tooltip(child=t("labels.track_length"),
                                        input_id='param_length_gen_track',
                                        tooltip_id="faq_param_length_gen_track",
                                        tooltip_text=t("tooltips.track_length"),
                                        input_type='number',
                                        input_value=t("default_input_value.track_length")),

                input_case_with_tooltip(child=t("labels.nb_tracks"),
                                        input_id='param_nb_tracks_gen_track',
                                        tooltip_id="faq_param_nb_tracks_gen_track",
                                        tooltip_text=t("tooltips.nb_tracks"),
                                        input_type='number',
                                        input_value=t("default_input_value.nb_tracks")),

                input_case(child=t("labels.filename"),
                           input_id='param_filename_gen_track',
                           input_type='text',
                           input_value=t("default_input_value.filename_data")),

            ], width=4),

            # Graph area
            dbc.Col(children=[
                dcc.Graph(id='profile-plot'),
            ], width=8),
        ]),

        html.Br(),
        html.Br(),

        # Button row
        dbc.Row([
            # Show plot profile button
            dbc.Col(children=[button(child=t("button.show_profile"),
                                     id_='show-profile-btn')],
                    width=2),

            # Generate track Button and Spinner Side by Side
            dbc.Col([
                html.Div([
                    dbc.Col(children=[
                        dcc.Store(id="start", data=""),
                        dcc.Store(id="complete", data=""),
                        dbc.Button(dbc.Spinner(
                            html.Span(children="Generate tracks",
                                      id="loading_generate")),
                            id="start-gen-tracks-btn"),
                    ], width=3),
                ]),
                html.Div(id='gen-tracks-output'),
                dcc.Download(id="download-csv2"),
            ]),
        ]),
    ]),
    )
