import numpy as np

from dash import  html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate

from ..analysis.analysis_track import recommend_acquisition

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
            print("Calculating acquisition parameters...")
            results = recommend_acquisition(k_est = params[2],
                               protein_length = params[0],
                               suntag_length = params[1],
                               nb_full_prot=params[3],
                               samples_per_ramp=params[4]
                               )
            
            return repr(results)
        raise PreventUpdate
    