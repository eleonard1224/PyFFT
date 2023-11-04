"""Create Test Functions"""

import numpy as np
from numpy.random import normal
from typing import List


def sinc_signal(**kwargs) -> np.ndarray:
    t = kwargs["time"]
    return np.sinc(t)


def hat_signal(
    tleft: float = -1.0, tright: float = 1.0, T: float = 0.25, N: int = 1024
) -> np.ndarray:
    t = np.linspace(tleft, tright, N)
    signal = np.zeros(N)
    unit_domain = np.where((t >= -T) & (t <= T))
    signal[unit_domain] = 1.0
    return signal
