"""Plotting Functions"""

import numpy as np
import matplotlib.pyplot as plt


# plotting functions
def plot_signal(time: np.ndarray, signal: np.ndarray) -> None:
    plt.figure(1)
    plt.plot(time, signal)


def plot_transform(frequency: np.ndarray, dft: np.ndarray, fft: np.ndarray) -> None:
    plt.figure(2)
    plt.plot(frequency, dft, color="limegreen", linewidth=10.0, label="DFT")
    plt.plot(frequency, fft, color="mediumvioletred", linewidth=2.5, label="FFT")
    plt.legend(loc="upper right")
