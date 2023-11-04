"""Calculate Fourier Transforms"""

import argparse
import time as timer
import numpy as np

from transforms import *
from signals import *
from plot import *


# parse command-line inputs
parser = argparse.ArgumentParser()
parser.add_argument(
    "--test-function",
    default="sinc-signal",
    help="Signal in time domain",
)
parser.add_argument(
    "--time-sample-endpoints", default="-5.0 5.0", help="Time sample endpoints"
)
parser.add_argument(
    "--number-of-time-samples", default="1024", help="Number of time samples"
)
args = parser.parse_args()

# read in command-line inputs
function = args.test_function.replace("-", "_")
time_sample_endpoints = list(
    map(lambda x: float(x), args.time_sample_endpoints.split(" "))
)
tleft = time_sample_endpoints[0]
tright = time_sample_endpoints[1]
N = int(args.number_of_time_samples)

# generate test function in time domain
time = np.linspace(tleft, tright, N)
dt = time[1] - time[0]
signal = eval(function)(time=time)

# plot signal in time domain
plot_signal(time, signal)

# get frequencies
f = get_frequencies(dt, N)

# calculate discrete fourier transform in O(N^2) time
t1 = timer.time()
dft = DFT(signal, N, dt)
t2 = timer.time()

# calculate fast fourier transform in O(N log N) time
t3 = timer.time()
fft = FFT(signal, N, dt)
t4 = timer.time()

# print timings
print(f"DFT Time = {round(t2-t1, 3)} sec, FFT Time = {round(t4-t3, 3)} sec")

# plot frequency signal
plot_transform(f, np.abs(dft), np.abs(fft))

# show plots
plt.show()
