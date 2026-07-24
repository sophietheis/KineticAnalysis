from dash import Input, Output, State
from dash.exceptions import PreventUpdate

import plotly.graph_objs as go

from .utils import empty_error_figure
from ..analysis.contribution import calculate_contribution
from ..plots.plots import fig_contribution

def register_callbacks(app):
    @app.callback(
        Output('equation-plot', 'figure'),
        Input('show-contribution-btn', 'n_clicks'),
        State('param_prot_length', 'value'),  #0
        State('param_suntag_length', 'value'),  #1
        State('param_nb_suntag', 'value'),   #2
        State('param_elongation_rate', 'value'),  #3
        State('param_initiation_rate', 'value'),   #4
        State('param_tau', 'value'),  #5
    )
    def update_plot(n_clicks, *params):
        """
        This function generate and plot an example for the simulation.
        """
        if n_clicks:
            try:
                L = float(params[0])
                s = float(params[2])
                N = int(float(params[1])/s)
                k = float(params[3])
                c = float(params[4])
                tau = float(params[5])

                print(L, s, N, k, c, tau)

                # Calculate contribution
                (term1,
                 term2,
                 term3,
                 term4,
                 term5) = calculate_contribution(L, N, k, c, s, tau)

                figure = fig_contribution(tau, term1, term2, term3, term4,
                                          term5)

                return figure

            except Exception as e:
                print(e)
                return empty_error_figure()
        raise PreventUpdate
