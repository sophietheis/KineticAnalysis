import multipletau
import numpy as np
import pandas as pd

from scipy import optimize


def read_csv_file(f):
    """
    Read csv file of trajectories.
    Drop first lines.
    Switch some type columns (turn it into numeric values)
    """

    datas = pd.read_csv(f)
    datas.drop(index=[0, 1, 2], inplace=True)
    datas['FRAME'] = pd.to_numeric(datas["FRAME"])
    datas['POSITION_X'] = pd.to_numeric(datas["POSITION_X"])
    datas['POSITION_Y'] = pd.to_numeric(datas["POSITION_Y"])
    datas['TRACK_ID'] = pd.to_numeric(datas["TRACK_ID"])
    datas['MEAN_INTENSITY_CH1'] = pd.to_numeric(datas["MEAN_INTENSITY_CH1"])
    datas['POSITION_T'] = pd.to_numeric(datas["POSITION_T"])
    datas.drop("MANUAL_SPOT_COLOR", axis=1, inplace=True)
    datas = datas.dropna(axis=0)

    return datas


def fit_function(x, t, c):
    """
    function : (((T-x)/(c*T**2)) * np.heaviside((T-x),0)))
    x : intensity signal
    T : residence time T = M/k where M = protein aa size and k is the elongation rate
    c : translation initiation rate
    """
    return ((t - x) / (c * t ** 2)) * np.heaviside((t - x), 0)


def autocorrelation(y, delta_t=0.5, normalize=True):
    """
    Perform auto correlation

    Parameters
    ----------
    y : intensity signal
    delta_t : time between two images
    normalize : default True

    """
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


def fit_autocorrelation(x, y, func_=fit_function, method='lm', protein_size=1500):
    """
    Fit autocorrelation curve with func_
    Parameters
    ----------
    x, y: x and y values of autocorrelation curve
    func_  : function to fit
    method : method of fit resolution
    protein_size: in aa in order to calculation the elongation rate
    """

    popt, pcov = optimize.curve_fit(
        func_,
        x[1:],
        y[1:],
        method=method)

    elongation_r = protein_size / popt[0]
    translation_init_r = popt[1]

    return elongation_r, translation_init_r


def single_track_analysis(datas, 
                          id_track=0, 
                          delta_t=0.5, 
                          protein_size=1500, 
                          normalise_intensity=2 ** 16 * 100):
    x = (datas[datas.TRACK_ID == id_track].sort_values('FRAME')['POSITION_T'].values -
         min(datas[datas.TRACK_ID == id_track].sort_values('FRAME')['POSITION_T'].values))
    y = (datas[datas.TRACK_ID == id_track].sort_values('FRAME')['MEAN_INTENSITY_CH1'].values / normalise_intensity)

    x_auto, y_auto = autocorrelation(y, delta_t, True)

    elongation_r, translation_init_r = fit_autocorrelation(x_auto, y_auto, fit_function, protein_size=protein_size)

    return x, y, x_auto, y_auto, elongation_r, translation_init_r
