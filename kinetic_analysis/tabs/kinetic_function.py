import multipletau
import warnings

import numpy as np
import pandas as pd

from scipy import optimize
import scipy.signal
import scipy.io.wavfile


def autocorrelation(y, delta_t=0.5, normalize=True, mm=None):
    """
    Perform auto correlation.

    Parameters
    ----------
    y : list
        intensity signal
    delta_t : float
        time between two images
    normalize : bool
        normalize the result to the square of the average input signal and
        the factor M-k; default value : True
    mm : int
        defines the number of points on one level, must be an even integer
    """
    if mm is None:
        mm = int(len(y) / 2 - 1)
    if (mm % 2) == 0:
        autocor = multipletau.autocorrelate(
            y,
            m=mm,
            deltat=delta_t,
            normalize=normalize)
    else:
        autocor = multipletau.autocorrelate(
            y,
            m=mm + 1,
            deltat=delta_t,
            normalize=normalize)

    return autocor.flatten()[0::2], autocor.flatten()[1::2]


def fit_function(x, t, c):
    """
    Function used in the autocorrelation fit

    Parameters
    ----------
    x : float
        intensity signal
    t : float
        residence time T = M/k where M = protein aa size and k is the
    elongation rate
    c : float
        translation initiation rate

    Returns
    -------

    Description
    -----------
    The function used is function : (((T-x)/(c*T**2)) * np.heaviside((T-x),0))

    """
    return ((t - x) / (c * t ** 2)) * np.heaviside((t - x), 0)


import sympy as sp


def fit_function_string(equation):
    # Define symbols
    x, t, c = sp.symbols("x t c")

    # Sympify using sympy's full namespace
    expr = sp.sympify(equation, locals=sp.__dict__)

    # Convert the sympy expression to a callable function
    func_ = sp.lambdify((x, t, c), expr, modules=["numpy"])

    return func_


# old fit_autocorrelation
def fit_autocorrelation_original(x, y, func_=fit_function, method='lm',
                                 protein_size=1500, first_dot=True):
    """
    Fit autocorrelation curve with func_
    Parameters
    ----------
    x, y: x and y values of autocorrelation curve
    func_  : function to fit
    method : method of fit resolution
    protein_size: in aa in order to calculation the elongation rate
    first_dot : bool, take the account the first dot in the analysis
    """
    if not first_dot:
        x = x[1:]
        y = y[1:]
    # print("original method")
    popt, pcov = optimize.curve_fit(func_,
                                    x,
                                    y,
                                    method=method)

    elongation_r = protein_size / popt[0]
    translation_init_r = 1 / popt[1]

    return elongation_r, translation_init_r, np.sqrt(np.diag(pcov))


# old fit_autocorrelation_v2
def fit_autocorrelation_linear(x, y, protein_size=1200):
    """
    Fit autocorrelation using a linear method.

    Parameters
    ----------
    x : list
        time value
    y : list
        intensity fluorescence value
    protein_size : int
        size of the protein in amino acid

    Returns
    -------
    elongation_r : float
    translation_init_r : float

    Description
    -----------
    Fit a linear equation in the first part of the curve.
    First step is to find a sign change in the curve to extract the first
    part of the curve.
    Use ax+b equation fit on the beginning of the curve.
    """
    # Find the position in the list t where the sign change
    ysign = np.sign(np.array(np.diff(y)))
    signchange = ((np.roll(ysign, 1) - ysign) != 0).astype(int)
    signchange[0] = 0
    if len(np.where(signchange == 1)[0]) == 0:
        t = -1
    else:
        t = np.where(signchange == 1)[0][0]

    elongation_r = protein_size / x[t]
    if len(x[:t]) < 2:
        return -1, -1, [-1, -1]
    res_fit = np.polyfit(x[:t], y[:t], 1)
    translation_init_r = (res_fit[1] * x[t])
    return elongation_r, translation_init_r, [-1, -1]


def single_track_analysis(df,
                          id_track=0,
                          delta_t=0.5,
                          protein_size=1500,
                          normalise_intensity=1,
                          normalise_auto=True,
                          mm=None,
                          lowpass_=False,
                          cutoff=50,
                          rtol=1e-4,
                          method="original",
                          force_analysis=False,
                          first_dot=True,
                          simulation=False):
    """
    Analysis of one track inside a dataframe.

    Parameters
    ----------
    df : pd.df
        dataframe that contains tracks
    id_track : int
        id of the track that will be analysed
    delta_t : float
        time between two time point in sec
    protein_size: int
        size of the protein (+ suntag) in amino acid
    normalise_intensity : float
        value used for normalised the intensity, default value : 1
    normalise_auto : bool
        normalise for the autocorrelation.  normalise the result to the
        square of the average input signal and the factor M-k, default
        value : True
    mm : int
        defines the number of points on one level, must be an even integer,
        default value : None
    lowpass_ : bool
        low pass filter the data, default value : False
    cutoff : int
        frequency for the low pass filter, default value : 50
    rtol : float
        to check if time is continuous , default value : 1e-4
    method : str
        choose the method of the analysis, "linear" or "original"
    force_analysis : bool
        force the analysis even if criteria not reach, default value : False
    first_dot : bool
        use the first point in the analysis, default value :True
    simulation : bool
        define if the track come from a simulation (True) or an experiment
        (False), default value : False

    Returns
    -------
    x : np.array
        list of time point of the track
    y : np.array
        list of fluorescent intensity of the track
    x_auto : np.array
        list of time point of the autocorrelation
    y_auto : np.array
        list of G(t) of the autocorrelation
    elongation_r : float
        estimate elongation rate
    translation_init_r : float
        estimate translation_rate
    perr : float
        estimate error

    Description
    -----------
    This function analyse one track to extract the estimate elongation and
    initiation rates.
    Columns names needs to be:
    - "TRACK_ID", to extract one track with its ID
    - "FRAME", list of int corresponding to the i time point
    - "MEAN_INTENSITY_CH1", correspond to the fluorescence intensity
    if the dataframe doesn't have these names, use rename_columns function
    to rename column(s).
    """

    # Extract time point and multiply by delta_t to get the real time of
    # each frame
    x = (df[df["TRACK_ID"] == id_track].sort_values('FRAME')[
             'FRAME'].values -
         min(df[df["TRACK_ID"] == id_track].sort_values('FRAME')[
                 'FRAME'].values))
    if not simulation:
        x = x * delta_t
    # Extract intensity value
    y = (df[df["TRACK_ID"] == id_track].sort_values('FRAME')[
             'MEAN_INTENSITY_CH1'].values / normalise_intensity)

    # Check if time is continuous and fix it if gap not too big
    if not check_continuous_time(x, delta_t, rtol=rtol):
        # print("Time not continuous")
        times_diff = np.diff(x)[
            np.where(np.isclose(np.diff(x), delta_t, rtol=rtol) == False)]
        if (times_diff < (5 * delta_t)).all():
            # fix the time difference if it misses less than 5 points
            # print("to fix")
            i = 0
            while i < (len(x) - 1):
                if np.round(x[i] - x[i + 1], decimals=2) > delta_t:
                    x = x[:i + 1] + [(x[i] + x[i + 1]) / 2] + x[i + 1:]
                i += 1
        else:
            # print("not fix")
            if not force_analysis:
                return np.repeat(np.nan, 7)
            else:
                warnings.warn("Analysis is forced for track " + str(id_track),
                              UserWarning)

    # Apply a low pass filter
    if lowpass_:
        filtered = lowpass(y, cutoff, 1 / delta_t)
        y = filtered

    # Perform the autocorrelation
    x_auto, y_auto = autocorrelation(y, delta_t, normalise_auto, mm)

    # Apply the method of analysis
    if method == "original":
        (elongation_r,
         translation_init_r,
         perr) = fit_autocorrelation_original(x_auto,
                                              y_auto,
                                              fit_function,
                                              protein_size=protein_size,
                                              first_dot=first_dot)
    elif method == "linear":
        (elongation_r,
         translation_init_r,
         perr) = fit_autocorrelation_linear(x_auto,
                                            y_auto,
                                            protein_size=protein_size)
    else:
        (elongation_r, translation_init_r, perr) = np.nan, np.nan, np.nan

    return x, y, x_auto, y_auto, elongation_r, translation_init_r, perr


def check_continuous_time(x, dt, rtol=0.001):
    """
    Check if the time points of a track is continuous, i.e. there is no
    missing point in the tracking

    Parameters
    ----------
    x : list
        time points list
    dt : float
        time expected between two points
    rtol : float
        tolerance

    Returns
    -------
    boolean, True if track is continuous, else False
    """
    return np.allclose(np.diff(x), dt, rtol=rtol)


def lowpass(data: np.ndarray, cutoff: float, sample_rate: float,
            poles: int = 5):
    sos = scipy.signal.butter(poles, cutoff, 'lowpass', fs=sample_rate,
                              output='sos')
    filtered_data = scipy.signal.sosfiltfilt(sos, data)
    return filtered_data


def calculate_msd(x, y, z):
    """
    Calculate mean square displacement of the track.

    Parameters
    ----------
    x : np.array
        position in x_axis
    y : np.array
        position in y_axis
    z : np.array
        position in z_axis

    Returns
    -------
    msd : float
        mean square displacement value
    """
    r = (x ** 2 + y ** 2 + z ** 2) ** 0.5
    diff = np.diff(r)
    diff_sq = diff ** 2
    msd = [np.mean(diff_sq[0:i]) for i in range(1, len(diff_sq))]
    return msd
