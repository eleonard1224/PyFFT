"""Calculate Fourier Transforms"""

from __future__ import annotations
import sys
import numpy as np
import itertools
from typing import List


class Complex:
    """Complex Number Class"""

    def __init__(self, real: float = 0.0, imaginary: float = 0.0):
        self.real = real
        self.imaginary = imaginary
        self.r = np.sqrt(self.real * self.real + self.imaginary * self.imaginary)
        self.theta = (
            0.0 if real == 0 and imaginary == 0 else np.arccos(self.real / self.r)
        )

    def __add__(self, c: Complex) -> Complex:
        real = self.real + c.real
        imaginary = self.imaginary + c.imaginary
        return Complex(real, imaginary)

    def __sub__(self, c: Complex) -> Complex:
        real = self.real - c.real
        imaginary = self.imaginary - c.imaginary
        return Complex(real, imaginary)

    def __mul__(self, c: Complex) -> Complex:
        real = self.real * c.real - self.imaginary * c.imaginary
        imaginary = self.imaginary * c.real + self.real * c.imaginary
        return Complex(real, imaginary)

    def __pow__(self, c: float) -> Complex:
        real = np.power(self.r, c) * np.cos(c * self.theta)
        imaginary = np.power(self.r, c) * np.sin(c * self.theta)
        return Complex(real, imaginary)

    @staticmethod
    def exp(a: float, b: float) -> Complex:
        """Calculates exp(a+ib)"""
        real = np.exp(a) * np.cos(b)
        imaginary = np.exp(a) * np.sin(b)
        return Complex(real, imaginary)

    @staticmethod
    def convert_complex_list_to_numpy(F: List[Complex]) -> np.ndarray:
        N = len(F)
        f = np.zeros(N, dtype=np.complex_)
        for i, FF in enumerate(F):
            f[i] = complex(FF.real, FF.imaginary)
        return f


def get_frequencies(dt: float, N: int) -> np.ndarray:
    f = np.zeros(N)
    f[:] = np.arange(N) / (N * dt)
    return f - (1 / (2 * dt))


def DFT(f: np.ndarray, N: int, dt: float) -> np.ndarray:
    """O(N^2) Discrete Fourier Transform Implementation"""
    if not np.log2(N).is_integer():
        sys.exit("Number of samples must be a power of two")
    fc = [Complex(f[k], 0.0) for k in range(N)]
    F = [Complex(0, 0) for n in range(N)]
    W = Complex.exp(0.0, -2.0 * np.pi / float(N))
    for n in range(N):
        m = (n + (N // 2)) % N  # return output in the right order
        for k in range(N):
            F[m] = F[m] + fc[k] * pow(W, float(n * k))
    F = [Complex(dt, 0) * F[i] for i in range(N)]
    return Complex.convert_complex_list_to_numpy(F)


def swap(f: np.ndarray, i: int, j: int) -> np.ndarray:
    tempi = f[i]
    tempi1 = f[i + 1]
    tempj = f[j]
    tempj1 = f[j + 1]
    f[i] = tempj
    f[i + 1] = tempj1
    f[j] = tempi
    f[j + 1] = tempi1
    return f


def bit_reversal(f: np.ndarray, N: int) -> np.ndarray:
    NN = 2 * N
    j = 0
    for i in np.arange(0, NN, 2):
        if j > i:
            f = swap(f, i, j)
        m = N
        while m >= 2 and j > (m - 1):
            j = j - m
            m = m // 2
        j = j + m
    return f


def FFT(f: np.ndarray, N: int, dt: float) -> np.ndarray:
    """O(N log N) Fast Fourier Transform Implementation"""
    if not np.log2(N).is_integer():
        sys.exit("Number of samples must be a power of two")
    f = np.array(list(itertools.chain(*[[f[i], 0.0] for i in range(N)])))
    f = bit_reversal(f, N)
    NN = 2 * N
    n = 2
    while n < NN:
        theta = 2 * np.pi / n
        alpha = 2 * np.sin(theta / 2) * np.sin(theta / 2)
        beta = np.sin(theta)
        wr = 1.0
        wi = 0.0
        m = 2 * n
        for k in np.arange(0, n, 2):
            for i in np.arange(k, NN, m):
                j = i + n
                tempr = wr * f[j] - wi * f[j + 1]
                tempi = wr * f[j + 1] + wi * f[j]
                f[j] = f[i] - tempr
                f[j + 1] = f[i + 1] - tempi
                f[i] = f[i] + tempr
                f[i + 1] = f[i + 1] + tempi
            wrtemp = wr - (alpha * wr + beta * wi)
            witemp = wi - (alpha * wi - beta * wr)
            wr = wrtemp
            wi = witemp
        n = 2 * n

    # arrange f for plotting
    for i in np.arange(0, N, 2):
        j = i + N
        f = swap(f, i, j)

    # convert f into standard numpy complex form
    f = np.array([complex(f[i], f[i + 1]) for i in np.arange(0, NN, 2)])
    return dt * f
