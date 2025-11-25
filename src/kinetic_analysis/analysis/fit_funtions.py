import numpy as np
from math import factorial

def fit_function_exact(x, k, c, N, M):
    denominator = (c / k) ** 2 * ((N * (N + 1)) / 2 + N * M) ** 2

    term1 = c / k * np.exp(-k * x) * sum(
        (N - n) * (N - n + 1) * (((2 * N) + n + 1) / 6) * (
                    ((k * x) ** n) / np.float128(factorial(n))) for n in
        range(0, int(N))) / denominator

    term2 = c / k * np.exp(-k * x) * N * sum(
        n * ((2 * N - n + 1) / 2) * ((k * x) ** n) / (
            np.float128(factorial(int(n)))) for n in
        range(0, int(N))) / denominator
    term3 = c / k * np.exp(-k * x) * (N ** 2) * ((N + 1) / 2) * sum(
        ((k * x) ** n) / (np.float128(factorial(int(n)))) for n in
        range(int(N), int(M))) / denominator
    term4 = c / k * np.exp(-k * x) * N * sum(n * ((1 + n) / 2) * (
                ((k * x) ** (M + N - n)) / (
            np.float128(factorial(int(M + N - n))))) for n in
                                             range(1, int(N))) / denominator

    term5 = c / k * np.exp(-k * x) * (N ** 2) * sum(
        (M - n) * (((k * x) ** n) / (np.float128(factorial(int(n))))) for n in
        range(0, int(M))) / denominator

    out = term1 + term2 + term3 + term4 + term5

    return out

def fit_function_epitope(x, k, c, N):
    return (k / c) * (2 / 3) * (1 / ((N * (N + 1)) ** 2)) * np.exp(-k * x) * (
        sum((N - n) * (N - n + 1) * ((2 * N) + n + 1) * (
                    ((k * x) ** n) / (np.float128(factorial(n)))) for n in
            range(0, int(N)))
    )

def fit_function_approx(x, t, c):
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
