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

![signal](https://github.com/eleonard1224/PyFFT/assets/45740519/bc4e65ea-947b-4d3a-9cec-2d2e03fb2b4b) 
![transform](https://github.com/eleonard1224/PyFFT/assets/45740519/4b194ddc-8e47-4ac1-85d5-c98c832c07e1)

## Timings

| N        | DFT (sec) | FFT (sec) |
| -------- | --------- | --------- |
| 256      | 0.75      | 0.003     |
| 512      | 2.86      | 0.01      |
| 1024     | 11.68     | 0.01      |
| 2048     | 43.9      | 0.03      |
