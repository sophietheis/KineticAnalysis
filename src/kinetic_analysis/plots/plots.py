import numpy as np

import matplotlib.pyplot as plt

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from ..analysis.fit_functions import function_epitope, function_exact, function_approx


def fig_update_background(figure, n_rows=1, ncols=1, width=1000, height=800):
    for i in range(1, n_rows+1):
        for j in range(1, ncols+1):
            figure.update_xaxes(mirror=True,
                                ticks='outside',
                                showline=True,
                                linecolor='black',
                                gridcolor='lightgrey',
                                row=i,
                                col=j)
            figure.update_yaxes(mirror=True,
                                ticks='outside',
                                showline=True,
                                linecolor='black',
                                gridcolor='lightgrey',
                                row=i,
                                col=j)

    figure.update_layout(width=width,
                         height=height,
                         plot_bgcolor="white"
                         )
    return figure


def fig_analyse_track(x, y,
                      x_fix, y_fix,
                      x_auto, y_auto,
                      y_fit,  dt,
                      figure=None):
    if figure is None:
        figure = make_subplots(rows=3,
                               cols=1,
                               subplot_titles=["track profile",
                                               "autocorrelation",
                                               "residuals"
                                               ])

        # plot track profile
        figure.add_trace(go.Scatter(x=x * dt,
                                    y=y,
                                    mode="lines",
                                    name="Intensity",
                                    line_color="#1A51DB"),  # Blue
                         row=1,
                         col=1,
                         )

        figure.add_trace(go.Scatter(x=x_fix * dt,
                                    y=y_fix,
                                    mode="lines+markers",
                                    name="Corrected Intensity",
                                    line_color="#DBA41A"),  # Orange
                         row=1,
                         col=1,
                         )
        figure.update_xaxes(title_text='Time (sec)', row=1, col=1)
        figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

        # plot autocorrelation
        figure.add_trace(go.Scatter(x=x_auto,
                                    y=y_auto,
                                    mode="lines+markers",
                                    name="Autocorrelation",
                                    line_color="#000000"),  # Black
                         row=2,
                         col=1),

        figure.add_trace(
            go.Scatter(x=x_auto,
                       y=y_fit,
                       mode="lines+markers",
                       name="Fit",
                       line_color="#B80909"),  # Red
            row=2,
            col=1),

        figure.update_xaxes(title_text='Time delay (tau)',
                            row=2,
                            col=1)
        figure.update_yaxes(title_text='G(tau)', row=2, col=1)

        # plot residuals
        figure.add_trace(go.Scatter(x=x_auto,
                                    y=np.repeat(0, len(x_auto)),
                                    mode="lines",
                                    name="",
                                    line_color="#aeb6c2"),  # grey
                         row=3,
                         col=1)

        figure.add_trace(go.Scatter(x=x_auto,
                                    y=y_auto - y_fit[:len(y_auto)],
                                    mode="markers",
                                    name="residuals",
                                    line_color="#000000"),  # black
                         row=3,
                         col=1)

        figure.update_xaxes(title_text='X data',
                            row=3,
                            col=1)
        figure.update_yaxes(title_text='Delta', row=3, col=1)

        fig_update_background(figure, n_rows=3, ncols=1 )

    return figure


def fig_contribution(tau, term1, term2, term3, term4, term5, figure=None):
    if figure is None:
        figure = make_subplots(rows=2,
                               cols=1,
                               subplot_titles=(
                                   'Autocorrelation function profile',
                                   'Autocorrelation function profile ('
                                   'percentage)',
                               ))

    sum_term = term1 + term2 + term3 + term4 + term5

    colors = ["#F2A5A2", "#A2C8F2", "#4891E5", "#1A63B7", "#F2CDA2", "#000000"]
    names = ["Stemloop term", "Crossterm1", "Crossterm2", "Crosstem3",
             "Post stemloop term", "Total"]
    terms = [term1, term2, term3, term4, term5, sum_term]
    line_type = np.repeat("solid", 5)
    line_type = np.append(line_type, "dash")

    for i in range(len(names)):
        # Plot autocorrelation profile
        figure.add_trace(go.Scatter(x=np.arange(tau),
                                    y=terms[i].astype(float),
                                    mode='lines',
                                    line_color=colors[i],
                                    line={'dash': line_type[i]},
                                    name=names[i],
                                    legendgroup=names[i],
                                    showlegend=True),
                         row=1,
                         col=1)

        # Plot autocorrelation profile percentage
        figure.add_trace(go.Scatter(x=np.arange(tau),
                                    y=(terms[i] * 100 / sum_term).astype(
                                        float),
                                    mode='lines',
                                    line_color=colors[i],
                                    line={'dash': line_type[i]},
                                    name=names[i],
                                    legendgroup=names[i],
                                    showlegend=False),
                         row=2,
                         col=1)

    figure.update_xaxes(matches='x2', row=1)
    figure.update_xaxes(matches='x2', row=2)
    

    figure.update_xaxes(title_text='Tau (sec)', row=1, col=1)
    figure.update_yaxes(title_text='G(tau)', row=1, col=1)

    figure.update_xaxes(title_text='Tau (sec)', row=2, col=1)
    figure.update_yaxes(title_text='G(tau)(%)', row=2, col=1)

    figure = fig_update_background(figure, n_rows=2, ncols=1 )

    figure.update_layout(legend_title_text='G(T) subterm')

    return figure

def fig_equation(M, N, k, c, tau, figure=None):
    colors = ["#31e2ec", "#a56712", "#769e71", "#d801b9", 
              "#e3369d", "#0e4c0d", "#57b87b", "#4292d0",
              "#ff8633", "#fb9fca", "#fdbf6f", "#e31a1c", 
                "#b2df8a", "#33a02c",]
    
    if figure is None:
        figure = make_subplots(rows=1,
                                cols=3,
                                subplot_titles=(
                                    'Epitope function profile',
                                    'Exact function profile',
                                    'Post stemloop function profile',
                                ))

    x_auto = np.arange(0, int(tau), 0.1)
    nb_curves = 10
    cpt = 0
    for i in range(1, int(tau), int(tau/nb_curves)):
        if i!=1:
            i=i-1
        i=int(i)
        y_fit = function_epitope(x_auto[::i], k, c, N)

        figure.add_trace(go.Scatter(x=x_auto[::i],
                                     y=y_fit,
                                     mode='lines',
                                     line_color=colors[cpt],
                                     name=i,
                                     legendgroup=i,
                                     showlegend=True),
                        row=1, col=1)
        
        y_fit = function_exact(x_auto[::i], k, c, N, M).astype(float)

        figure.add_trace(go.Scatter(x=x_auto[::i],
                                    y=y_fit,
                                    mode='lines',
                                    line_color=colors[cpt],
                                    name=i, 
                                    legendgroup=i,
                                    showlegend=False),
                         row=1, col=2)



        y_fit = function_approx(x_auto[::i],N/k, c)
        figure.add_trace(go.Scatter(x=x_auto[::i],
                                    y=y_fit,
                                    mode='lines',
                                    line_color=colors[cpt],
                                    name=i,
                                    legendgroup=i,
                                    showlegend=False),
                        row=1, col=3)
        cpt+=1

    figure.update_xaxes(matches='x2', col=1)
    figure.update_xaxes(matches='x2', col=2)
    figure.update_xaxes(matches='x2', col=3)
    figure.update_layout(legend_title_text='Time step (sec)')
    figure = fig_update_background(figure, n_rows=1, ncols=3, width=1200, height=600)
    return figure


def fig_generate_track(x_profile, y_profile, x_track, y_track,
                       y_number, figure=None):
    if figure is None:
        figure = make_subplots(rows=3,
                               cols=1,
                               subplot_titles=(
                                   'One protein fluo profile',
                                   'One track fluo profile',
                                   'Number of translation'))
    # Plot one protein profile
    figure.add_trace(go.Scatter(x=x_profile, y=y_profile,
                                mode='lines',
                                name='Profile one prot'),
                     row=1,
                     col=1)
    figure.update_xaxes(title_text='Time (sec)', row=1, col=1)
    figure.update_yaxes(title_text='Fluorescence', row=1, col=1)

    # Plot one track
    figure.add_trace(go.Scatter(x=x_track, y=y_track,
                                mode='lines',
                                name='Profile track'),
                     row=2,
                     col=1)
    figure.update_xaxes(title_text='Time (sec)', row=2, col=1)
    figure.update_yaxes(title_text='Fluorescence', row=2, col=1)

    # Plot number of translation
    figure.add_trace(go.Scatter(x=x_track, y=y_number,
                                mode='lines',
                                name='Number of protein being '
                                     'translated'),
                     row=3,
                     col=1)
    figure.update_xaxes(title_text='Time (sec)', row=3, col=1)
    figure.update_yaxes(title_text='Number of translation', row=3,
                        col=1)
    figure.update_xaxes(matches='x2', row=2)
    figure.update_xaxes(matches='x2', row=3)
    figure.update_xaxes(matches=None, row=1)
    
    figure = fig_update_background(figure, 3)
    return figure

def plot_ribosome(df_single_prot, df_polysome, prot_size, suntag_size, figure=None):
    """
    Plot the ribosome density based on the single protein and polysome fluorescence data.

    Parameters
    ----------
    df_single_prot : pd.DataFrame
        Dataframe containing the single protein fluorescence data.
    df_polysome : pd.DataFrame
        Dataframe containing the polysome fluorescence data.
    prot_size : int
        Length of the protein of interest in amino acids.
    suntag_size : int
        Length of the tag in amino acids.

    Returns
    -------
    figure : plotly.graph_objs.Figure
        Figure object containing the ribosome density plot.
    """
    # mean_intensity_single, df_polysome = calculate_ribosome_density(
    #     df_single_prot, df_polysome, "intensity", "intensity", prot_size, suntag_size)

    
    if figure is None:
        figure = make_subplots(rows=1,
                                cols=3,
                                subplot_titles=(
                                    'Single protein fluo distribution',
                                    'Polysome fluo distribution',
                                    'Ribosome density distribution'))

    figure.add_trace(go.Histogram(df_single_prot,
                                  nbins=50),
        row=1,
        col=1)
    
    figure.add_trace(go.Histogram(df_polysome,
                                      nbins=50),
            row=1,
            col=2)

    figure.update_layout(
        title='Ribosome Density Estimation',
        xaxis_title='Sample Index',
        yaxis_title='Ribosome Density (ribosomes per amino acid)',
        template='plotly_white'
    )

    return figure
