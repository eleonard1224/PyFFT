# PyFFT

This repository contains Python code for a class-based implementation of the $O(N^{2})$ Discrete Fourier Transform as well as code which implements the $O(N \log{N})$ Cooley-Tukey Fast Fourier Transform algorithm given that the number of input samples is a power of two.

## Environment Setup

```
conda create --name pyfft python=3.8
conda activate pyfft
cd [PyFFT-root]
pip install -r requirements.txt
```

## Sample Results

## Timings

| N        | DFT (sec) | FFT (sec) |
| -------- | --------- | --------- |
| 256      | 0.75      | 0.003     |
| 512      | 2.86      | 0.01      |
| 1024     | 11.68     | 0.01      |
| 2048     | 43.9      | 0.03      |