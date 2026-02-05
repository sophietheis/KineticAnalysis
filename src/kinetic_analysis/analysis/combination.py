import random
import numpy as np
import pandas as pd


def combine_n_tracks(list_track=[]):
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
