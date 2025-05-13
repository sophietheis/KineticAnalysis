import warnings

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
    datas = pd.read_csv(f,
                        sep=None,
                        engine="python",
                        index_col=0)
    return datas


def read_csv_file_v1(f):
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
