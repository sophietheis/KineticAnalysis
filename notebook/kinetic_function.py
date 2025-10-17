import multipletau
import warnings

import numpy as np
import pandas as pd

from scipy import optimize
import scipy.signal
import scipy.io.wavfile


def generate_profile(prot_length,
                     suntag_length,
                     nb_suntag,
                     fluo_one_suntag,
                     translation_rate,
                     retention_time=0,
                     suntag_pos="begin",
                     step=0.1,
                     noise=False,
                     noise_std=1):
    """
    Generate fluorescence profile of one protein.

    Parameters
    ----------
    prot_length : int
        length of the protein in amino acid
    suntag_length : int
        length of the suntag in amino acid
    nb_suntag : int
        number of suntag repetition
    fluo_one_suntag : int
        fluorescence intensity of one suntag
    translation_rate : int
        translation rate of the protein in aa/sec
    retention_time : float, default 0
        length of time the protein remains on the translation site in sec
    suntag_pos : str, "begin" or "end", default "begin"
        position of the suntag, before of after the protein
    step : float, default 0.1
        time step between two point in sec
    noise : bool, default False
        add noise to the signal
    noise_std : float, default 1
        std of the normal distribution

    Returns
    -------
    x : list of time point
    y : list of fluorescent intensity

    Description
    -----------
    This function generate the fluorescent profile of one protein.
    The profile can have retention time, change the suntag position (before or
    after the protein).
    Time to "create" one protein is :
    $ (prot_length + suntag_length)/translation_rate + retention_time
    The assumption made for the fluorescence profile is there is no
    fluorescence intensity change during protein creation. During suntag
    creation, fluorescence intensity increase linearly according to
    the number of suntag and the translation rate.
    Examples
    --------
    """

    if suntag_pos == "begin":
        suntag_pos = 0
    elif suntag_pos == "end":
        suntag_pos = -1
    else:
        raise ValueError("suntag_pos value can only be \"begin\" or \"end\"")

    prot_tot_length = prot_length + suntag_length

    x = np.arange(prot_tot_length / translation_rate + retention_time,
                  step=step, )
    y = (nb_suntag * fluo_one_suntag) / (
            suntag_length / translation_rate) * np.arange(
        suntag_length / translation_rate,
        step=0.1)

    if suntag_pos == 0:
        y_prim = np.repeat(y[-1], len(x) - len(y))
        y = np.concatenate([y, y_prim])
    elif suntag_pos == -1:
        y_prim = np.repeat(0, len(x) - len(y))
        y = np.concatenate([y_prim, y])
    if noise:
        n = np.random.normal(0, noise_std, len(x))
        y = y + n

    return x, y


def generate_track(prot_length,
                   suntag_length,
                   nb_suntag,
                   fluo_one_suntag,
                   translation_rate,
                   binding_rate,
                   retention_time=0,
                   suntag_pos="begin",
                   step=0.1,
                   length=6000):
    """
    Generate track according to one protein translation dynamics

    Parameters
    ----------
    prot_length : int
        length of the protein in amino acid
    suntag_length : int
        length of the suntag in amino acid
    nb_suntag : int
        number of suntag repetition
    fluo_one_suntag : int
        fluorescence intensity of one suntag
    translation_rate : int
        translation rate of the protein in aa/sec
    binding_rate : float
        probability to start a new protein in 1/sec
    retention_time : float, default 0
        length of time the protein remains on the translation site in sec
    suntag_pos : str, "begin" or "end", default "begin"
        position of the suntag, before of after the protein
    step : float, default 0.1
        time step between two point in sec
    length : int
        length of the track in sec

    Returns
    -------
    x_global : list of time point
    y_global : list of fluorescent intensity
    y_start_prot : list of number of protein in translation

    Description
    -----------
    This function generate one track of a translation site to mimic in vivo
    translation profile.
    It is based on one protein translation profile.

    """
    x, y = generate_profile(prot_length,
                            suntag_length,
                            nb_suntag,
                            fluo_one_suntag,
                            translation_rate,
                            retention_time,
                            suntag_pos,
                            step)

    # global signal
    x_global = np.arange(length, step=step)
    y_global = np.zeros(len(x_global))
    y_start_prot = np.zeros(len(x_global))

    n_rand = np.random.rand(len(x_global))
    for i in range(len(x_global)):
        # random number between 0 and 1
        if n_rand[i] < (binding_rate * step):
            if i > (len(x_global) - len(x)):
                y_global[i:i + len(x)] += y[:len(y_global[i:i + len(x)])]
                y_start_prot[i:i + len(x)] += 1
            else:
                y_global[i:i + len(x)] += y
                y_start_prot[i:i + len(x)] += 1

    # Remove the first time points
    x_global = x_global[2000:] - 200  # -200 to start time at 0
    y_global = y_global[2000:]
    y_start_prot = y_start_prot[2000:]

    return x_global, y_global, y_start_prot


def read_csv_file(f):
    """
    Read csv file of trajectories.
    Drop first lines.
    Switch some type columns (turn it into numeric values)
    """

    datas = pd.read_csv(f, sep=None, engine="python")
    # datas.drop(index=[0, 1, 2], inplace=True)
    datas['FRAME'] = pd.to_numeric(datas["FRAME"])
    datas['POSITION_X'] = pd.to_numeric(datas["POSITION_X"])
    datas['POSITION_Y'] = pd.to_numeric(datas["POSITION_Y"])
    try:
        datas['TRACK_ID'] = pd.to_numeric(datas["TRACK_ID"])
    except:
        pass
    datas['MEAN_INTENSITY_CH1'] = pd.to_numeric(datas["MEAN_INTENSITY_CH1"])
    datas['POSITION_T'] = pd.to_numeric(datas["POSITION_T"])
    datas.drop("MANUAL_SPOT_COLOR", axis=1, inplace=True)
    datas = datas.dropna(axis=0)

    return datas


def read_csv_file_v2(f):
    """
    Read csv file of trajectories.
    Drop first lines.
    Switch some type columns (turn it into numeric values)
    """

    datas = pd.read_csv(f, sep=None, engine="python")
    datas.columns = datas.iloc[0]
    datas.columns = [x.upper() for x in datas.columns]
    datas.drop([0, 1, 2], axis=0, inplace=True)
    datas.reset_index(inplace=True)
    datas.drop(index=[0, 1, 2], inplace=True)
    datas['FRAME'] = pd.to_numeric(datas["FRAME"])
    datas['X'] = pd.to_numeric(datas["X"])
    datas['Y'] = pd.to_numeric(datas["Y"])
    datas['Z'] = pd.to_numeric(datas["Z"])
    datas['TRACK ID'] = pd.to_numeric(datas["TRACK ID"])
    datas['MEAN INTENSITY CH1'] = pd.to_numeric(datas["MEAN INTENSITY CH1"])
    datas['T'] = pd.to_numeric(datas["T"])
    datas.drop("MANUAL SPOT COLOR", axis=1, inplace=True)
    datas = datas.dropna(axis=0)

    return datas


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


def rename_columns(df, old_columns, new_columns):
    """
    Rename columns name

    Parameters
    ----------
    df : dataframe
    old_columns : list
        list of columns name to be replaced
    new_columns : list
        list of columns new name used to be replaced
    """
    if len(old_columns) != len(new_columns):
        raise "lengths of old_columns is different of new_columns"
    for i in range(len(old_columns)):
        df.rename(columns={old_columns[i]: new_columns[i]}, inplace=True)


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


def test(abc):
    return abc
