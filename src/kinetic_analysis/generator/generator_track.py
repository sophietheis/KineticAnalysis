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
                       footprint=0,
                       retention_time=0,
                       suntag_pos="begin",
                       noise=False,
                       noise_std=0,
                       step=0.1,
                       length=6000,
                       remove_point_beginning=2000
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
    footprint : int
        ribosome exclusion size in amino acids (default None, no exclusion)
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
    remove_point_beginning : int
        how many point are removed at the beginning of the track

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
    if footprint == -1:
        # Generate one protein profile
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

        # Generate the track based on the protein profile and the binding rate
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

    else:
        prot_tot_length = prot_length + suntag_length
        x_global = np.arange(length, step=step)
        y_global = np.zeros(len(x_global))
        y_start_prot = np.zeros(len(x_global))

        positions = []          # aa position, ordered leading (idx 0) -> trailing (idx -1)
        retention_timers = []   # None while elongating, float while retained/terminating

        n_rand = np.random.rand(len(x_global))

        for i in range(len(x_global)):
            # Elongation, leading ribosome first, blocked if ribosome ahead is < footprint away
            for idx in range(len(positions)):
                if retention_timers[idx] is not None:
                    continue
                ahead = positions[idx - 1] if idx > 0 else np.inf
                if ahead - positions[idx] > footprint:
                    positions[idx] += translation_rate * step
                    if positions[idx] >= prot_tot_length:
                        positions[idx] = prot_tot_length
                        retention_timers[idx] = retention_time if retention_time > 0 else 0.0

            # Retention countdown / termination
            keep = []
            for idx in range(len(positions)):
                if retention_timers[idx] is not None:
                    retention_timers[idx] -= step
                    if retention_timers[idx] > 0:
                        keep.append(idx)
                else:
                    keep.append(idx)
            positions = [positions[j] for j in keep]
            retention_timers = [retention_timers[j] for j in keep]

            # Initiation, blocked if first `footprint` aa occupied
            can_initiate = (len(positions) == 0) or (positions[-1] > footprint)
            if can_initiate and n_rand[i] < (binding_rate * step):
                positions.append(0.0)
                retention_timers.append(None)

            # Record fluorescence and active ribosome count (for density estimation)
            total_fluo = 0.0
            for pos in positions:
                if pos < suntag_length:
                    total_fluo += (nb_suntag * fluo_one_suntag) * (pos / suntag_length)
                else:
                    total_fluo += nb_suntag * fluo_one_suntag
            y_global[i] = total_fluo
            y_start_prot[i] = len(positions)


    # Remove the first time points
    x_global = (x_global[remove_point_beginning:] -
                (remove_point_beginning * step)) 
    y_global = y_global[remove_point_beginning:]
    y_start_prot = y_start_prot[remove_point_beginning:]


    return x_global, y_global, y_start_prot


def generate_tracks(n,
                    prot_length,
                    suntag_length,
                    nb_suntag,
                    fluo_one_suntag,
                    translation_rate,
                    binding_rate,
                    footprint=-1,
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
    footprint : int
        ribosome exclusion size in amino acids (default None, no exclusion)
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
                                                              footprint,
                                                              retention_time,
                                                              suntag_pos,
                                                              noise,
                                                              noise_std,
                                                              step,
                                                              length)

        frame = pd.DataFrame({"FRAME": x_global,
                            "MEAN_INTENSITY_CH1": y_global,
                            "TRACK_ID": i,
                            "RETENTION_TIME": retention_time,
                            "N_RIBOSOME": y_start_prot
                                  })
        datas = frame if first_time else pd.concat([datas, frame], ignore_index=True)
        first_time = False
        
    return datas
