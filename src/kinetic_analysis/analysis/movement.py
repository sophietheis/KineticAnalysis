import numpy as np


def msd_calculation(x, y, z):
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
    # diff = np.sqrt(np.diff(x)**2 + np.diff(y)**2 + np.diff(z)**2)
    # diff = np.concatenate([[0], diff])

    r = (x ** 2 + y ** 2 + z ** 2) ** 0.5
    diff = np.diff(r)
    diff = np.concatenate([[0], diff])
    diff_sq = diff ** 2
    msd = [np.mean(diff_sq[0:i]) for i in range(1, len(diff_sq)+1)]

    return msd


def rsd_calculation(x, y, z):

    data = np.array([x, y, z]).T
    # Compute squared displacement for each time step
    sd = [(i[0] - x[0]) ** 2 + (i[1] - y[0]) ** 2 + (i[2] - z[0]) ** 2 for i in
          data]
    # Compute the cumulative average of SD to get MSD at each time step
    rsd = np.cumsum(sd) / np.arange(1, len(sd) + 1)
    return rsd
