from uuid import uuid4

from dash import dcc, Input, Output, State
from dash.exceptions import PreventUpdate

import plotly.graph_objs as go

from ..generator.generator_track import (generate_one_track,
                                         generate_tracks,
                                         generate_profile)
from ..plots.plots import fig_generate_track

from ..tabs.app_function import (browse_directory)


def register_callbacks(app):
    @app.callback(
        Output("directory-output", "children"),
        [Input("select_directory", "n_clicks")],
    )
    def select_directory(n_clicks):
        """
        Select the directory in which to save the file.
        """
        return browse_directory(n_clicks, 'directory_generation', app)

    @app.callback(
        Output('profile-plot', 'figure'),
        Input('show-profile-btn', 'n_clicks'),
        State('param_prot_length', 'value'),
        State('param_suntag_length', 'value'),
        State('param_nb_suntag', 'value'),
        State('param_fluo_one_suntag', 'value'),
        State('param_translation_rate', 'value'),
        State('param_initiation_rate', 'value'),
        State('param_retention_time', 'value'),
        State('param_pos_suntag', 'value'),
        State('param_noise', 'value'),
        State('param_dt', 'value'),
        State('param_length', 'value'),
        State('param_footprint', 'value'),
    )
    def update_profile_plot(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """
        if n_clicks:
            try:
                # Generate profile
                noise = False
                if params[8] > 0:
                    noise = True

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
                print(e)
                return {
                    'data': [],
                    'layout': go.Layout(title='Error', xaxis={'title': 'Time'},
                                        yaxis={'title': 'Fluorescence'})
                }
        raise PreventUpdate

    # make the button works
    @app.callback(
        Output("start-gen-tracks-btn", "disabled", allow_duplicate=True),
        Output("start", "data"),
        Input("start-gen-tracks-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def register_start(n):
        if n:
            return True, str(uuid4())
        raise PreventUpdate

    @app.callback(
        Output("start-gen-tracks-btn", "disabled", allow_duplicate=True),
        Input("complete", "data"),
        State("start", "data"),
        prevent_initial_call=True,
    )
    def enable_button(complete_value, start_value):
        return complete_value != start_value

    @app.callback(
        Output('gen-tracks-output', 'children'),
        Output('download-csv2', 'data'),
        Output("loading_generate", "children"),
        Output("complete", "data"),
        Input('start-gen-tracks-btn', 'n_clicks'),
        Input("start", "data"),
        State('param_prot_length', 'value'),
        State('param_suntag_length', 'value'),
        State('param_nb_suntag', 'value'),
        State('param_fluo_one_suntag', 'value'),
        State('param_translation_rate', 'value'),
        State('param_initiation_rate', 'value'),
        State('param_retention_time', 'value'),
        State('param_pos_suntag', 'value'),
        State('param_noise', 'value'),
        State('param_dt', 'value'),
        State('param_length', 'value'),
        State('param_nb_tracks', 'value'),
        State('param_filename', 'value'),
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
                return (f"Error: {str(e)}",
                        None,
                        "Generate tracks",
                        data)
        raise PreventUpdate
