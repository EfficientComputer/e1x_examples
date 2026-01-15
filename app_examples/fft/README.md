# FFT

### Specification

* Implements a complex FFT with bit-reversed output.
* Input is an array 16-bit Q15 format integers. Each integer is a single PCM sample.
* Output is an array of 16-bit Q15 format integers. Each integer represents a bin's value.
* Because of limited precision, FFT size may not exceed 4096. For this example, it is set to 4096 bins.
* Element error compared to the expected ideal must be lower than 10 as defined by `abs(expected[i] - output[i])`. The ideal is found using the unmodified KissFFT library.

### FFT Overview

FFT, short for Fast Fourier Transform, converts between the time domain and the frequency domain
in DSP. FFTs are said to have n-bins, where n is a fixed number depending on the task. Each bin
corresponds to a slice of frequencies in the frequency domain or a PCM sample in the time domain.