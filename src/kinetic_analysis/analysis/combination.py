import random
import numpy as np
import pandas as pd

from .analysis_track import autocorrelation


def combine_n_tracks(list_track=[]):
    """
    Concatenates tracks with intensity offset matching. NOTE: this
    introduces spurious correlation at track junctions and biases fitted
    k downward as more tracks are combined (see Discussion). For kinetic
    estimation, prefer ensemble_autocorrelation / combine_tracks_df_ensemble.
    Kept here for direct comparison against the ensemble-averaged method.
    """
    if len(list_track) == 0:
        return None

    y_concat = np.empty(shape=[0])
    i = 0
    for lt in list_track:
        if len(lt) != 0:
            if i == 0:
                y_concat = np.concatenate([y_concat, lt])
                i += 1
            else:
                last_val = y_concat[-1]
                diff_ = last_val - lt[0]

                y_concat = np.concatenate([y_concat, (lt + diff_)[1:]])

    return y_concat


def reservoir_sampling(sampled_size, total_num, nb_combination):
    """

    :param sampled_size: size of the new subgroup
    :type sampled_size:
    :param total_num: size of the original group
    :type total_num:
    :param nb_combination: number of subgroup created
    :type nb_combination:
    :return:
    :rtype:
    """
    pools = []
    for k in range(nb_combination):
        pool = []
        for i in range(0, total_num):
            if i < sampled_size:
                pool.append(i)
            else:
                r = random.randint(0, i)
                if r < sampled_size:
                    pool[r] = i
        pools.append(pool)
    return pools


def combine_tracks_df(df, n_tracks, n_size):
    """

    :param df: dataframe containing all datas
    :type df:
    :param n_tracks: number of new tracks created
    :type n_tracks:
    :param n_size: number of tracks in new created tracks
    :type n_size:
    :return:
    :rtype:
    """
    new_df = pd.DataFrame(columns=["TRACK_ID", "MEAN_INTENSITY_CH1","FRAME",
                                   "ORIGINAL_IDs"])

    id_tracks = np.unique(df["TRACK_ID"])
    list_combinations = reservoir_sampling(n_size, len(id_tracks), n_tracks)

    a=0
    for lc in list_combinations:
        intensities = []
        real_id_tracks = []
        for l in lc :
            intensities.append(df[df["TRACK_ID"] == id_tracks[l]][
                "MEAN_INTENSITY_CH1"].to_numpy())
            real_id_tracks.append(id_tracks[l])
        combined_tracks = combine_n_tracks(intensities)

        new_df = pd.concat([new_df,
                            pd.DataFrame(
                                {"TRACK_ID": np.repeat(a, len(combined_tracks)),
                                 "FRAME": np.arange(len(combined_tracks)),
                                 'MEAN_INTENSITY_CH1': combined_tracks,
                                 'ORIGINAL_IDs': np.repeat([str(
                                     real_id_tracks)],
                                                           len(combined_tracks),
                                                           axis=0)},
                            )
                            ], ignore_index=True)
        a+=1
    return new_df


def ensemble_autocorrelation(list_y, delta_t=0.5, normalize=True, mm=None):
    """
    Ensemble-averaged G(tau) from independent tracks, avoiding the
    concatenation artifact of combine_n_tracks. Each track only
    contributes within-track pairs; results are pooled weighted by
    (T_i - tau).

    Parameters
    ----------
    list_y : list of np.ndarray
        List of intensity traces to compute autocorrelation for.
    delta_t : float, optional
        Time step between frames, by default 0.5.
    normalize : bool, optional
        Whether to normalize the autocorrelation, by default True.
    mm : float, optional
        Mean intensity to use for normalization, by default None.   

    Outputs
    -------
    all_tau : np.ndarray
        Array of time lags (tau) for the autocorrelation.
    g : np.ndarray
        Ensemble-averaged autocorrelation values corresponding to all_tau.
    
    """
    all_tau, num, denom = None, None, None
    for y in list_y:
        if len(y) < 4:
            continue
        tau_i, g_i = autocorrelation(y, delta_t=delta_t, normalize=normalize, mm=mm)
        weight = (len(y) - np.arange(len(tau_i))).astype(float)
        if all_tau is None:
            all_tau, num, denom = tau_i, g_i * weight, weight
        else:
            n = min(len(tau_i), len(all_tau))
            num, denom, all_tau = num[:n] + (g_i * weight)[:n], denom[:n] + weight[:n], all_tau[:n]
    return all_tau, num / denom


def combine_tracks_df_ensemble(df, n_tracks, n_size, delta_t=0.5, normalize=True, mm=None):
    """
    Ensemble-averaged alternative to combine_tracks_df. Returns tau/G(tau)
    pooled across n_size independently sampled tracks, ready to pass
    directly to fit_autocorrelation_* — no synthetic concatenated
    intensity trace is created.
    """
    id_tracks = np.unique(df["TRACK_ID"])
    list_combinations = reservoir_sampling(n_size, len(id_tracks), n_tracks)

    results = []
    for lc in list_combinations:
        intensities = [df[df["TRACK_ID"] == id_tracks[l]]["MEAN_INTENSITY_CH1"].to_numpy() for l in lc]
        tau, g = ensemble_autocorrelation(intensities, delta_t, normalize, mm)
        results.append({"tau": tau, "G": g, "ORIGINAL_IDs": [id_tracks[l] for l in lc]})
    return results