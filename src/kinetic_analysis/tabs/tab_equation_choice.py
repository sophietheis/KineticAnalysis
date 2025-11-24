import numpy as np

from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc


import plotly.graph_objs as go
from plotly.subplots import make_subplots


from ..analysis.contribution import calculate_contribution


def layout():
    return (html.Div([
        # Explanation at the begining of the page
        dbc.Row([
            html.P([
                "In this tab, you will be able to visualise the "
                "contribution of each term in the equation in order "
                "to help you choose the best method of analysis. ",
                html.Br(),
                "This is based on the work of Larson et al. ",
                html.Br(),
                ]),
            html.Br(),
        ]),
        # Define parameter of the simulation
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.P(["Protein length (aa) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_prot_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}
                           ),
                    dbc.Tooltip("Length of the protein in amino acid. This "
                                "value will be added to suntag length.",
                                target="faq_param_prot_length"),
                    dcc.Input(id='param_prot_length', type='number', value=490,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Suntag length (aa) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_suntag_length",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Length of the suntag in amino acid. This "
                                "value will be added to protein length.",
                                target="faq_param_suntag_length"),
                    dcc.Input(id='param_suntag_length', type='number', value=796,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Number of suntag ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_nb_suntag",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height": "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Number of suntag repetition",
                                target="faq_param_nb_suntag"),
                    dcc.Input(id='param_nb_suntag', type='number', value=32,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Estimate elongation rate (aa/sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_translation_rate",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip("Elongation rate.",
                                target="faq_param_translation_rate"),
                    dcc.Input(id='param_elongation_rate', type='number',
                              value=0.3,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["Estimate initiation rate (ribosome/sec) ",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_initiation_rate",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={
                               "height":
                                   "auto",
                               "margin-bottom": "auto"}),
                    dbc.Tooltip("Initiation rate.",
                                target="faq_param_initiation_rate"),
                    dcc.Input(id='param_initiation_rate', type='number',
                              value=0.015,
                              style={'width': '200px'}),
                ]),
                html.Div([
                    html.P(["tau",
                            html.Span(className="fas fa-question-circle",
                                      id="faq_param_dt",
                                      style={"cursor": "pointer",
                                             "marginLeft": "5px"})],
                           style={"height":
                                      "auto",
                                  "margin-bottom": "auto"}),
                    dbc.Tooltip(
                        "Maximum length",
                        target="faq_param_dt"),
                    dcc.Input(id='param_tau', type='number', value=150,
                              style={'width': '200px'}),
                ]),
                html.Br(),
            ], width=3),

            # Show plot contribution
            dbc.Col([
                dbc.Button('Show Contribution', id='show-contribution-btn',
                           className="mr-1"),
                dcc.Graph(id='equation-plot'),

            ], width=5)
        ]),
    ]),
    )


## Callbacks
def register_callbacks(app):
    @app.callback(
        Output('equation-plot', 'figure'),
        Input('show-contribution-btn', 'n_clicks'),
        State('param_prot_length', 'value'),  #0
        State('param_suntag_length', 'value'),  #1
        State('param_nb_suntag', 'value'),   #2
        State('param_elongation_rate', 'value'),  #3
        State('param_initiation_rate', 'value'),   #4
        State('param_tau', 'value'),  # 5
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

                print(L,s,N,k,c,tau)

                # Calculate contribution
                (term1,
                 term2,
                 term3,
                 term4,
                 term5) = calculate_contribution(L, N, k, c, s, tau)

                sum_term = term1 + term2 + term3 + term4 + term5

                # Create the figure
                figure = make_subplots(rows=2,
                                       cols=1,
                                       subplot_titles=(
                                           'Autocorrelation function profile',
                                           'Autocorrelation function profile (percentage)',
                                           ))

                # Plot autocorrelation profile
                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=term1.astype(float),
                                            mode='lines',
                                            line_color="#F2A5A2",
                                            name='stemloop term'),
                                 row=1,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=term2.astype(float),
                                            mode='lines',
                                            line_color="#A2C8F2",
                                            name='Crossterm1'),
                                 row=1,
                                 col=1)
                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=term3.astype(float),
                                            mode='lines',
                                            line_color="#4891E5",
                                            name='Crossterm2'),
                                 row=1,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=term4.astype(float),
                                            mode='lines',
                                            line_color="#1A63B7",
                                            name='Crossterm3'),
                                 row=1,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=term5.astype(float),
                                            mode='lines',
                                            line_color="#F2CDA2",
                                            name='post stem loop term'),
                                 row=1,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=(sum_term).astype(float),
                                            mode='lines',
                                            line_color="#000000",
                                            line={'dash':"dash"},
                                            name='Total'),
                                 row=1,
                                 col=1)

                # Plot autocorrelation profile percentage
                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=(term1*100/sum_term).astype(
                                                float),
                                            mode='lines',
                                            line_color="#F2A5A2",
                                            name='stemloop term'),
                                 row=2,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=(term2*100/sum_term).astype(
                                                float),
                                            mode='lines',
                                            line_color="#A2C8F2",
                                            name='Crossterm1'),
                                 row=2,
                                 col=1)
                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=(term3*100/sum_term).astype(
                                                float),
                                            mode='lines',
                                            line_color="#4891E5",
                                            name='Crossterm2'),
                                 row=2,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=(term4*100/sum_term).astype(
                                                float),
                                            mode='lines',
                                            line_color="#1A63B7",
                                            name='Crossterm3'),
                                 row=2,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=(term5*100/sum_term).astype(float),
                                            mode='lines',
                                            line_color="#F2CDA2",
                                            name='post stem loop term'),
                                 row=2,
                                 col=1)

                figure.add_trace(go.Scatter(x=np.arange(tau),
                                            y=(sum_term*100/sum_term).astype(
                                                float),
                                            mode='lines',
                                            line_color="#000000",
                                            line={'dash': "dash"},
                                            name='Total'),
                                 row=2,
                                 col=1)

                figure.update_xaxes(title_text='Tau (sec)', row=1, col=1)
                figure.update_yaxes(title_text='G(tau)', row=1, col=1)

                figure.update_xaxes(title_text='Tau (sec)', row=2, col=1)
                figure.update_yaxes(title_text='G(tau)(%)', row=2, col=1)

                for i in range(1, 3):
                    figure.update_xaxes(mirror=True,
                                        ticks='outside',
                                        showline=True,
                                        linecolor='black',
                                        gridcolor='lightgrey',
                                        row=i,
                                        col=1)
                    figure.update_yaxes(mirror=True,
                                        ticks='outside',
                                        showline=True,
                                        linecolor='black',
                                        gridcolor='lightgrey',
                                        row=i,
                                        col=1)

                figure.update_layout(width=1000,
                                     height=800,
                                     plot_bgcolor = "white"
                                     )
                return figure

            except Exception as e:
                print(e)
                return {
                    'data': [],
                    'layout': go.Layout(title='Error', xaxis={'title': 'Time'},
                                        yaxis={'title': 'Fluorescence'})
                }
        raise PreventUpdate
