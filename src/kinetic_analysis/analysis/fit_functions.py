import numpy as np
from math import factorial

def function_exact(x, k, c, N, M):
    denominator = (c / k) ** 2 * ((N * (N + 1)) / 2 + N * M) ** 2

    term1 = c / k * np.exp(-k * x) * sum(
        (N - n) * (N - n + 1) * (((2 * N) + n + 1) / 6) * (
                    ((k * x) ** n) / np.float128(factorial(n))) for n in
        range(0, int(N)))

    term2 = c / k * np.exp(-k * x) * N * sum(
        n * ((2 * N - n + 1) / 2) * ((k * x) ** n) / (
            np.float128(factorial(int(n)))) for n in
        range(0, int(N)))

    term3 = c / k * np.exp(-k * x) * (N ** 2) * ((N + 1) / 2) * sum(
        ((k * x) ** n) / (np.float128(factorial(int(n)))) for n in
        range(int(N), int(M)))

    term4 = c / k * np.exp(-k * x) * N * sum(n * ((1 + n) / 2) * (
                ((k * x) ** (M + N - n)) / (
            np.float128(factorial(int(M + N - n))))) for n in
                                             range(1, int(N)))

    term5 = c / k * np.exp(-k * x) * (N ** 2) * sum(
        (M - n) * (((k * x) ** n) / (np.float128(factorial(int(n))))) for n in
        range(0, int(M)))

    out = (term1 + term2 + term3 + term4 + term5) / denominator

    return np.longdouble(out)

def function_epitope(x, k, c, N):
    return (k / c) * (2 / 3) * (1 / ((N * (N + 1)) ** 2)) * np.exp(-k * x) * (
        sum((N - n) * (N - n + 1) * ((2 * N) + n + 1) * (
                    ((k * x) ** n) / (np.float128(factorial(n)))) for n in
            range(0, int(N)))
    )

def function_approx(x, t, c):
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



#### A supprimer ?
def fit_function_linear(x, y):
    """
    Fit autocorrelation using a linear method.

    Parameters
    ----------
    x : list
        time value
    y : list
        aucorrelation value
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
        t_sign = -1
    else:
        t_sign = np.where(signchange == 1)[0][0]


    # find when the curve cross the x_axis
    y_sign_value = np.sign(y)
    t_xaxis = np.where(y_sign_value==-1)[0][0]
    # If the sign change happen in a negative value
    if t_sign < t_xaxis:
        t = t_xaxis
    else:
        t = t_sign
    # print(t_sign, t)
    # elongation_r = protein_size / x[t]
    if len(x[:t]) < 2:
        return -1, -1, [-1, -1]
    res_fit = np.polyfit(x[:t], y[:t], 1)
    # translation_init_r = (res_fit[1] * x[t])
    return(res_fit[0], res_fit[1], [np.nan, np.nan])