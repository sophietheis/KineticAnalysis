import warnings

import numpy as np
import pandas as pd


def generate_profile(prot_length,
                     suntag_length,
                     nb_suntag,
                     fluo_one_suntag,
                     translation_rate,
                     retention_time=0,
                     suntag_pos="begin",
                     step=0.1,
                     noise=False,
                     noise_std=0):
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
        step=step)

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


def generate_one_track(prot_length,
                       suntag_length,
                       nb_suntag,
                       fluo_one_suntag,
                       translation_rate,
                       binding_rate,
                       retention_time=0,
                       suntag_pos="begin",
                       noise=False,
                       noise_std=0,
                       step=0.1,
                       length=6000,
                       ):
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
    noise : bool, default False
        add noise to the signal
    noise_std : float, default 1
        std of the normal distribution
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
                            step,
                            noise,
                            noise_std
                            )

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


def generate_tracks(n,
                    prot_length,
                    suntag_length,
                    nb_suntag,
                    fluo_one_suntag,
                    translation_rate,
                    binding_rate,
                    retention_time=0,
                    suntag_pos="begin",
                    noise=False,
                    noise_std=0,
                    step=0.1,
                    length=6000):
    """
    Generate n tracks according to one protein translation dynamics

    Parameters
    ----------
    n : int
        number of tracks
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
    noise : bool, default False
        add noise to the signal
    noise_std : float, default 1
        std of the normal distribution
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
    first_time = True

    for i in range(n):
        x_global, y_global, y_start_prot = generate_one_track(prot_length,
                                                              suntag_length,
                                                              nb_suntag,
                                                              fluo_one_suntag,
                                                              translation_rate,
                                                              binding_rate,
                                                              retention_time,
                                                              suntag_pos,
                                                              noise,
                                                              noise_std,
                                                              step,
                                                              length)
        if first_time:
            datas = pd.DataFrame({"FRAME": x_global,
                                  "MEAN_INTENSITY_CH1": y_global,
                                  "TRACK_ID": i,
                                  "RETENTION_TIME": retention_time,
                                  })
            first_time = False
        else:
            datas = pd.concat([datas,
                               pd.DataFrame({"FRAME": x_global,
                                             "MEAN_INTENSITY_CH1": y_global,
                                             "TRACK_ID": i,
                                             "RETENTION_TIME": retention_time
                                             })],
                              ignore_index=True)

    return datas