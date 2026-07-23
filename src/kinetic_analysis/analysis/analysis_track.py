import multipletau
import lmfit

import numpy as np

from scipy import optimize

from .fit_functions import (function_exact,
                            function_approx,
                            function_epitope)

from .analyse_density import estimate_density


def correct_elongation_rate(k_fit, rho_bar):
    """
    Correct the elongation rate for ribosome queuing using the mean occupancy

    Parameters
    ----------
    k_fit : float
        fitted elongation rate
    rho_bar : float
        mean occupancy of the mRNA  

    Returns
    -------
    k_true : float
        corrected elongation rate

    """
    if rho_bar is None or rho_bar >= 1:
        return np.nan
    return k_fit / (1 - rho_bar)

def correct_initiation_rate(c_fit, rho_bar):
    """
    Correct the initiation rate for ribosome queuing using the mean occupancy

    Parameters
    ----------
    c_fit : float
        fitted initiation rate
    rho_bar : float
        mean occupancy of the mRNA  

    Returns
    -------
    c_true : float
        corrected initiation rate   

    """
    if rho_bar is None or rho_bar >= 1:
        return np.nan
    return c_fit / (1 - rho_bar)

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


def fit_autocorrelation_exact(x, y, M=56, N=32):
    model = lmfit.Model(function_exact)
    params = lmfit.Parameters()
    params.add('N', value=N, vary=False)
    params.add('M', value=M, vary=False)
    params.add('k', value=np.float128(0.6), min=1e-15)
    params.add('c', value=np.float128(0.1), min=1e-15)

    try:
        result = model.fit(y, params, x=x, nan_policy='raise')
    except ValueError as e:
        print("Sorry, due to the complexity of the equation, there is a "
              "chance that M or N are too big, and the number are too large "
              "to be processed now (cause of factorial).")
        print(e)
        return np.nan, np.nan, [np.nan, np.nan]

    return (result.params["k"].value, result.params["c"].value,
            [result.params["k"].stderr, result.params["c"].stderr])


def fit_autocorrelation_epitope(x, y, N=32):
    model = lmfit.Model(function_epitope)
    params = lmfit.Parameters()
    params.add('N', value=N, vary=False)
    params.add('k', value=np.float128(0.6), min=0)
    params.add('c', value=np.float128(0.1), min=0)

    result = model.fit(y, params, x=x, nan_policy='raise')

    return (result.params["k"].value, result.params["c"].value,
            [result.params["k"].stderr, result.params["c"].stderr])


def fit_autocorrelation_approx(x, y, method='lm'):
    """
    Fit autocorrelation curve with func_
    Parameters
    ----------
    x, y: x and y values of autocorrelation curve
    method : method of fit resolution
    """

    popt, pcov = optimize.curve_fit(function_approx,
                                    x,
                                    y,
                                    method=method)

    # elongation_r = protein_size / popt[0]
    # translation_init_r = 1 / popt[1]
    # return elongation_r, translation_init_r, np.sqrt(np.diag(pcov))
    return popt[0], popt[1], np.sqrt(np.diag(pcov))


def single_track_analysis(x,
                          y,
                          delta_t=0.5,
                          protein_size=1500,
                          suntag_size=800,
                          repetition_suntag=32,
                          normalise_auto=True,
                          mm=None,
                          method="exact",
                          simulation=False,
                          mean_n_ribosome=None,  
                          footprint=10,          
                          correct_queuing=False, ):
    """
    Analysis of one track inside a dataframe.

    Parameters
    ----------
    x : pd.df
        dataframe that contains tracks
    y : int
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
    rtol : float
        to check if time is continuous , default value : 1e-4
    method : str
        choose the method of the analysis, "linear" or "original"
    force_analysis : bool
        force the analysis even if criteria not reach, default value : False
    simulation : bool
        define if the track come from a simulation (True) or an experiment
        (False), default value : False
    func_ : function
        contains the function that is used for the fit

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
    :param suntag_size:
    :type suntag_size:
    """

    if not simulation:
        x = x * delta_t

    # Perform the autocorrelation
    x_auto, y_auto = autocorrelation(y, delta_t, normalise_auto, mm)

    one_suntag_size = (int(suntag_size/repetition_suntag))
    M = int(protein_size/one_suntag_size)
    N = repetition_suntag
    print(M, N)
    # Apply the method of analysis
    if method == "exact":
        (k, c, perr) = fit_autocorrelation_exact(x_auto,
                                                 y_auto,
                                                 N=N,
                                                 M=M)
        elongation_r = k*(suntag_size/repetition_suntag)
        translation_init_r = c

    elif method == "approx":
        (k, c, perr) = fit_autocorrelation_approx(x_auto,
                                                  y_auto,)

        elongation_r = M/k*one_suntag_size
        # translation_init_r = (1 / (c * k))
        translation_init_r = c
    elif method == "epitope":
        (k, c, perr) = fit_autocorrelation_epitope(x_auto,
                                                   y_auto,
                                                   N=N)

        elongation_r = k*one_suntag_size
        # translation_init_r = (1 / (c * k))
        translation_init_r = c
    else:
        print("No method choose")
        (k, c, elongation_r, translation_init_r, perr) = (np.nan, np.nan,
                                                          np.nan, np.nan,
                                                          [np.nan, np.nan])
    if correct_queuing:
        if mean_n_ribosome is None:
            raise ValueError("mean_n_ribosome must be provided when correct_queuing=True")
        rho_bar = estimate_density(mean_n_ribosome, protein_size, suntag_size, footprint)
        elongation_r = correct_elongation_rate(elongation_r, rho_bar)
        translation_init_r = correct_initiation_rate(translation_init_r, rho_bar)

    return x_auto, y_auto, k, c, elongation_r, translation_init_r, perr


def check_track_validity(df,
                         id_track,
                         delta_t,
                         normalise_intensity=1,
                         simulation=False,
                         rtol=1e-4,
                         nb_missing_point=5):
    # Extract time point and multiply by delta_t to get the real time of
    # each frame
    x = (df[df["TRACK_ID"] == id_track].sort_values('FRAME')[
             'FRAME'].values -
         min(df[df["TRACK_ID"] == id_track].sort_values('FRAME')[
                 'FRAME'].values))

    # Extract intensity value
    y = (df[df["TRACK_ID"] == id_track].sort_values('FRAME')[
             'MEAN_INTENSITY_CH1'].values / normalise_intensity)
    x_orig = x
    y_orig = y

    if not simulation:
        x = x * delta_t

    # Check if time is continuous and fix it if gap not too big
    if not check_continuous_time(x, delta_t, rtol=rtol):

        if (np.diff(x) < (nb_missing_point * delta_t)).all():
            # fix the time difference if it misses less than nb_missing_point
            i = 0
            while i < (len(x) - 1):

                if np.round(x[i + 1] - x[i], decimals=2) > delta_t:

                    # How many points to add
                    x_add = [x[i]+delta_t]
                    y_add = [(y[i]+y[i+1]) / 2]

                    while x_add[-1] < (x[i+1]-delta_t):
                        x_add.append(x_add[-1]+delta_t)
                        y_add.append((y_add[-1]+y[i+1])/2)

                    x = np.concatenate((x[:i+1],
                                        np.array(x_add),
                                        x[i+1:]))
                    y = np.concatenate((y[:i+1],
                                        np.array(y_add),
                                        y[i+1:]))
                i += 1
            i += 1
        else:
            print("Gap is too big - not fix")
            return False, x_orig, y_orig, x/delta_t, y

    return True, x_orig, y_orig, x/delta_t, y


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


def recommend_acquisition(k_est, protein_length, suntag_length, nb_full_prot=7, samples_per_ramp=24):
    """
    Recommend acquisition parameters based on estimated elongation rate and protein length.

    Parameters
    ----------
    k_est : float
        Estimated elongation rate in amino acids per second.
    protein_length : int
        Length of the protein in amino acids.
    suntag_length : int
        Length of the suntag in amino acids.
    nb_full_prot : int, optional
        Number of proteins produced entirely one after the other
    samples_per_ramp : int, optional
        Number of samples to take during the ramp phase (default is 24). 
    
    Returns
    -------
    dict
        A dictionary containing recommended acquisition parameters:
        - 'tau_c': Total translation time (protein + suntag) / k_est
        - 'tau_ramp': Translation time for the protein only / k_est
        - 'T_recommended': Total acquisition time to capture nb_full_prot times
        - 'dt_recommended': Recommended time interval between samples
        - 'dt_nyquist_limit': Nyquist limit for time interval 
    """
    tau_c = (protein_length + suntag_length) / k_est
    tau_ramp = suntag_length / k_est
    T_recommended = nb_full_prot * tau_c
    dt_recommended = tau_ramp / samples_per_ramp
    return {"tau_c": tau_c, "tau_ramp": tau_ramp,
            "T_recommended": T_recommended,
            "dt_recommended": dt_recommended,
            "dt_nyquist_limit": tau_ramp / 2,
            "n_points": int(T_recommended / dt_recommended)}