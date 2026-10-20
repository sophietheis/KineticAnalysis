import logging

from uuid import uuid4

from dash import dcc, Input, Output, State
from dash.exceptions import PreventUpdate

from .utils import empty_error_figure
from ..generator.generator_track import (generate_one_track,
                                         generate_tracks,
                                         generate_profile)
from ..plots.plots import fig_generate_track

logger = logging.getLogger(__name__)


def register_callbacks(app):

    @app.callback(
        Output('gentrack_plot_profile', 'figure'),
        Input('gentrack_button_profile', 'n_clicks'),
        State('gentrack_prot_length', 'value'),
        State('gentrack_epitope_length', 'value'),
        State('gentrack_epitope_repeat', 'value'),
        State('gentrack_epitope_fluo', 'value'),
        State('gentrack_elongation_rate', 'value'),
        State('gentrack_initiation_rate', 'value'),
        State('gentrack_retention_time', 'value'),
        State('gentrack_epitope_position', 'value'),
        State('gentrack_noise', 'value'),
        State('gentrack_dt', 'value'),
        State('gentrack_track_length', 'value'),
        State('gentrack_footprint', 'value'),
    )
    def update_profile_plot(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """
        if n_clicks:
            try:
                # Generate profile
                noise = params[8] > 0

                if params[11] == -1:
                    x_profile, y_profile = generate_profile(prot_length=int(params[0]),
                                                            suntag_length=int(params[1]),
                                                            nb_suntag=int(params[2]),
                                                            fluo_one_suntag=int(params[3]),
                                                            translation_rate=int(params[4]),
                                                            retention_time=int(params[6]),
                                                            suntag_pos=params[7],
                                                            noise=noise,
                                                            noise_std=float(params[8]),
                                                            step=float(params[9]))
                else:
                    x_profile, y_profile = None, None

                # Generate one track
                x, y, y_number = generate_one_track(prot_length=int(params[0]),
                                                    suntag_length=int(params[1]),
                                                    nb_suntag=int(params[2]),
                                                    fluo_one_suntag=int(params[3]),
                                                    translation_rate=int(params[4]),
                                                    binding_rate=float(params[5]),
                                                    footprint=int(params[11]),
                                                    retention_time=float(params[6]),
                                                    suntag_pos=params[7],
                                                    noise=noise,
                                                    noise_std=float(params[8]),
                                                    step=float(params[9]),
                                                    length=int(params[10])
                                                    )

                figure = fig_generate_track(x_profile, y_profile, x, y,
                                            y_number)

                return figure

            except Exception as e:
                logger.exception("Failed to generate profile/track plot", e)
                return empty_error_figure()
        raise PreventUpdate

    # make the button works
    @app.callback(
        Output("gentrack_button_start", "disabled", allow_duplicate=True),
        Output("gentrack_start", "data"),
        Input("gentrack_button_start", "n_clicks"),
        prevent_initial_call=True,
    )
    def register_start(n):
        if n:
            return True, str(uuid4())
        raise PreventUpdate

    @app.callback(
        Output("gentrack_button_start", "disabled", allow_duplicate=True),
        Input("gentrack_complete", "data"),
        State("gentrack_start", "data"),
        prevent_initial_call=True,
    )
    def enable_button(complete_value, start_value):
        return complete_value != start_value

    @app.callback(
        Output('gentrack_output', 'children'),
        Output('gentrack_download_csv', 'data'),
        Output("gentrack_loading", "children"),
        Output("gentrack_complete", "data"),
        Input('gentrack_button_start', 'n_clicks'),
        Input("gentrack_start", "data"),
        State('gentrack_prot_length', 'value'),
        State('gentrack_epitope_length', 'value'),
        State('gentrack_epitope_repeat', 'value'),
        State('gentrack_epitope_fluo', 'value'),
        State('gentrack_elongation_rate', 'value'),
        State('gentrack_initiation_rate', 'value'),
        State('gentrack_retention_time', 'value'),
        State('gentrack_epitope_position', 'value'),
        State('gentrack_noise', 'value'),
        State('gentrack_dt', 'value'),
        State('gentrack_track_length', 'value'),
        State('gentrack_nb_track', 'value'),
        State('gentrack_filename', 'value'),
        # background=True,
        # progress=[Ouput("progress-generate", "children")]
    )
    def start_generate_tracks(n_clicks, data,  *params):
        # Generate all tracks and save it
        if n_clicks:
            try:
                noise = False
                if params[8] > 0:
                    noise = True
                datas = generate_tracks(n=int(params[11]),
                                        prot_length=int(params[0]),
                                        suntag_length=int(params[1]),
                                        nb_suntag=int(params[2]),
                                        fluo_one_suntag=int(params[3]),
                                        translation_rate=int(params[4]),
                                        binding_rate=float(params[5]),
                                        retention_time=float(params[6]),
                                        suntag_pos=params[7],
                                        noise=noise,
                                        noise_std=float(params[8]),
                                        step=float(params[9]),
                                        length=int(params[10]),
                                        )

                output_path = str(params[12]) + ".csv"
                datas.to_csv(output_path, index=False)

                return ("Tracks generated and saved successfully!",
                        dcc.send_file(output_path),
                        "Generate tracks",
                        data)
            except Exception as e:
                logger.exception("Failed to generate tracks")
                return (f"Error: {str(e)}",
                        None,
                        "Generate tracks",
                        data)
        raise PreventUpdate
