import multipletau
import numpy as np
import pandas as pd

from scipy import optimize
import scipy.signal
import scipy.io.wavfile


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


def autocorrelation(y, delta_t=0.5, normalize=True, mm=None):
    """
    Perform auto correlation

    Parameters
    ----------
    y : intensity signal
    delta_t : time between two images
    normalize : default True

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
    print("original method")
    popt, pcov = optimize.curve_fit(
        func_,
        x,
        y,
        method=method)

    elongation_r = protein_size / popt[0]
    translation_init_r = popt[1]

    return elongation_r, translation_init_r

def fit_autocorrelation_v2(x, y, protein_size=1200):
    # fit ax+b equation at the begining, until curve reach 0
    # Find t
    print("linear method")
    ysign = np.sign(np.array(y))
    signchange = ((np.roll(ysign, 1) - ysign) != 0).astype(int)
    signchange[0] = 0
    t = np.where(signchange==1)[0][0]

    elongation_r = protein_size/x[t]
    if len(x[:t])<2:
        return -1, -1
    res_fit = np.polyfit(x[:t], y[:t], 1)
    print(res_fit)
    translation_init_r = 1/(res_fit[1]*x[t])

    return elongation_r, translation_init_r



def lowpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 5):
    sos = scipy.signal.butter(poles, cutoff, 'lowpass', fs=sample_rate, output='sos')
    filtered_data = scipy.signal.sosfiltfilt(sos, data)
    return filtered_data


def single_track_analysis(datas,
                          id_track=0,
                          delta_t=0.5,
                          protein_size=1500,
                          normalise_intensity=2 ** 16 * 100,
                          normalize_auto=True,
                          mm=None,
                          lowpass_=False,
                          cutoff=50,
                          rtol=1e-4,
                          method="original",
                          force_analysis=False):
    x = (datas[datas.TRACK_ID == id_track].sort_values('FRAME')['POSITION_T'].values -
         min(datas[datas.TRACK_ID == id_track].sort_values('FRAME')['POSITION_T'].values))
    y = (datas[datas.TRACK_ID == id_track].sort_values('FRAME')['MEAN_INTENSITY_CH1'].values / normalise_intensity)

    if not check_continuous_time(x, delta_t, rtol=rtol):
        times_diff = np.diff(x)[np.where(np.isclose(np.diff(x), delta_t, rtol=rtol) == False)]
        if (times_diff < 3 * delta_t).all():
            print("to fix")
            i=0
            while i < (len(x)-1):
                if np.round(x[i]-x[i+1], decimals=2)>delta_t:
                    x = x[:i+1] + [(x[i]+x[i+1])/2] + x[i+1:]
                i+=1
        else:
            print("not fix")
            if not force_analysis:
                return np.repeat(np.nan, 6)
            else:
                print("force analysis")

    if lowpass_:
        filtered = lowpass(y, cutoff, 1 / delta_t)
        y = filtered

    # y = y - np.mean(y)
    x_auto, y_auto = autocorrelation(y, delta_t, normalize_auto, mm)

    if method == "original":
        elongation_r, translation_init_r = fit_autocorrelation(x_auto, y_auto, fit_function, protein_size=protein_size)
    elif method == "linear":
        elongation_r, translation_init_r = fit_autocorrelation_v2(x_auto, y_auto, protein_size=protein_size)

    return x, y, x_auto, y_auto, elongation_r, translation_init_r


def check_continuous_time(x, dt, rtol=0.001):
    return np.allclose(np.diff(x), dt, rtol=rtol)
