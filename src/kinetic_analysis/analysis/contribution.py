import numpy as np
from math import factorial
from decimal import Decimal


def calculate_contribution(L, N, k, c, s, t):
    '''

    :param L:
    :type L:
    :param N:
    :type N:
    :param k:
    :type k:
    :param c:
    :type c:
    :return:
    :rtype:
    '''
    print("calculate_contribution")
    M = int(L / s)
    tau = np.arange(0, np.int32(t))

    denominator = (c / k) ** 2 * ((N * (N + 1)) / 2 + N * M) ** 2

    k = np.float128(k)
    c = np.float128(c)

    term1 = c / k * np.exp(-k * tau) * sum(
        (N - n) * (N - n + 1) * (((2 * N) + n + 1) / 6) * (
                    ((k * tau) ** n) / np.float128(factorial(n))) for n in
        range(0, int(N))) / denominator

    term2 = c / k * np.exp(-k * tau) * N * sum(
        n * ((2 * N - n + 1) / 2) * ((k * tau) ** n) / (
            np.float128(factorial(int(n)))) for n in
        range(0, int(N))) / denominator
    term3 = c / k * np.exp(-k * tau) * (N ** 2) * ((N + 1) / 2) * sum(
        ((k * tau) ** n) / (np.float128(factorial(int(n)))) for n in
        range(N, int(M))) / denominator
    term4 = c / k * np.exp(-k * tau) * N * sum(n * ((1 + n) / 2) * (
                ((k * tau) ** (M + N - n)) / (
            np.float128(factorial(int(M + N - n))))) for n in
                                               range(1, int(N))) / denominator

    term5 = c / k * np.exp(-k * tau) * (N ** 2) * sum(
        (M - n) * (((k * tau) ** n) / (np.float128(factorial(int(n))))) for n
        in range(0, int(M))) / denominator

    return term1, term2, term3, term4, term5