import logging 

import numpy as np

from dash import  html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate

from ..analysis.analysis_track import recommend_acquisition

logger = logging.getLogger(__name__)

## Callbacks
def register_callbacks(app):

    @app.callback(
        Output('output_acquisition', 'children'),
        Input('btn_calculate_acquisition', 'n_clicks'),
        State('param_prot_length_acquisition', 'value'),
        State('param_suntag_length_acquisition', 'value'),
        State('param_translation_rate_acquisition', 'value'),
        State('param_full_translated_protein_acquisition', 'value'),
        State('param_stem_loops_acquisition', 'value'),
    )
    def calculate_acquisition(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """
        if n_clicks:
            logger.debug("Calculating acquisition parameters...")
            try:
                results = recommend_acquisition(k_est = params[2],
                                protein_length = params[0],
                                suntag_length = params[1],
                                nb_full_prot=params[3],
                                samples_per_ramp=params[4]
                                )
            except (TypeError, ZeroDivisionError) as e:
                logger.exception("Failed to compute recommended acquisition parameters")
                return f"Could not compute acquisition parameters: {e}. Please check all fields are filled in."            
            new_line = '<br>'
            out = (f"Total time to translate (SunTag + protein) = {results['tau_c']:.2f} sec, \n" \
                  f"Time to translate SunTag part = {results['tau_ramp']:.2f} sec, \n" \
                  f"Total acquisition time = {results['T_recommended']:.2f} sec, \n"  \
                  f"Recommended frame rate = {results['dt_recommended']:.2f} sec, \n" \
                #   f"Nyquist limit = {results['dt_nyquist_limit']:.2f} sec, \n" \
                  f"Number of points = {int(results['T_recommended'] / results['dt_recommended']):.0f}, " 
            )

            logger.debug("Acquisition parameters calculated: %s", out)
            return html.Pre(out)
        raise PreventUpdate
    