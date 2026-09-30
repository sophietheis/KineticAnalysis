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

def plot_ribosome(df_single_prot, df_polysome, result, figure=None):
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
        figure =  make_subplots(rows=1,
                                cols=2,
                                # shared_xaxes=True,
                                # shared_yaxes=True,
                                subplot_titles=(
                                    'Fluorescent distribution',
                                    'Ribosome density per polysome',)
                                )

    figure.add_trace(go.Histogram(x=df_single_prot, nbinsx=50, name='Single protein'),
                     row=1, col=1, )
    
    figure.add_trace(go.Histogram(x=df_polysome, nbinsx=50, name='Polysome'),
                     row=1, col=1, )
    
    figure.update_xaxes(title_text='Fluorescence intensity', row=1, col=1)
    figure.update_yaxes(title_text='Count', row=1, col=1)

    figure.add_trace(go.Scatter(x=result["INTENSITY"],
                                y=result["ribosome_density"],
                                mode='markers',
                                name='Ribosome density'),
                     row=1, col=2)

    figure.update_xaxes(title_text='Polysome fluorescence intensity', row=1, col=2)
    figure.update_yaxes(title_text='Ribosome density (ribosomes/aa)', row=1, col=2)

    return figure

def plot_ribosome_density(prot_size, suntag_size, rib_density,fig=None):
    """
    Plot the ribosome density based on the single protein and polysome fluorescence data.

    Parameters
    ----------
    prot_size : int
        Length of the protein of interest in amino acids.
    suntag_size : int
        Length of the tag in amino acids.
    rib_density : float
        The calculated ribosome density.
    

    Returns
    -------
    figure : plotly.graph_objs.Figure
        Figure object containing the ribosome density plot.
    """

    if fig is None:
        fig = go.Figure()

    total = prot_size + suntag_size  
    rib_number = int(rib_density * total)  # Calculate the number of ribosomes based on density and total length
    step = total / rib_number 


    # --- Rectangles ---
    fig.add_shape(
        type="rect", x0=0, y0=0, x1=suntag_size, y1=2,
        fillcolor="red", opacity=0.5, line_width=0
    )
    fig.add_shape(
        type="rect", x0=suntag_size, y0=0, x1=suntag_size + prot_size, y1=2,
        fillcolor="blue", opacity=0.5, line_width=0
    )

    fig.add_annotation(x=0,    y=1, text="suntag", showarrow=False,
                    xanchor="left", font=dict(size=12))
    fig.add_annotation(x=suntag_size, y=1, text="prot",   showarrow=False,
                    xanchor="left", font=dict(size=12))

    # --- Circles (evenly spaced) ---

    fig.add_trace(go.Scatter(
        x=[i * step for i in range(rib_number)],
        y=[2.5] * rib_number,
        mode="markers",
        marker=dict(symbol="circle", size=12, color="black"),
        showlegend=False,
        hoverinfo="skip"
    ))

    # --- Ribosome bracket ---
    fig.add_trace(go.Scatter(
        x=[0, 0, step, step],
        y=[3.2, 3.6, 3.6, 3.2],
        mode="lines",
        line=dict(color="black", width=2),
        showlegend=False,
        hoverinfo="skip"
    ))
    fig.add_annotation(x=0, y=3.8, text="ribosome", showarrow=False,
                    xanchor="left", font=dict(size=12))

    # --- Footprint bracket ---
    tilt = 5
    fig.add_trace(go.Scatter(
        x=[tilt * step / 2, tilt * step / 2,
        (tilt + 2) * step / 2, (tilt + 2) * step / 2],
        y=[3.7, 4.1, 4.1, 3.7],
        mode="lines",
        line=dict(color="black", width=2),
        showlegend=False,
        hoverinfo="skip"
    ))
    fig.add_annotation(x=tilt * step / 2, y=4.2, text="footprint", showarrow=False,
                    xanchor="left", font=dict(size=12))

    # --- Layout ---
    fig.update_layout(
        xaxis=dict(range=[-100, total+100], title="Amino Acids", showgrid=False, zeroline=False),
        yaxis=dict(range=[-1, 6], visible=False),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=20, r=20, t=20, b=40),
        width=800,
        height=400
    )

    return fig

