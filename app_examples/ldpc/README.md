# LDPC

## Specification

* Performs LDPC encoding and decoding. Decoding is much more complex, and is the step for which the profile regions are captured. The code is a (336, 672) code copied from the 802.11 WiFi standard. The decoder uses min-sum approximation for belief propagation, and the provided test case is expected to pass after 2 iterations.

* Input (to decoder) is an int32 array of the log likelihood ratios for each bit in the codeword. An LLR of -inf means the bit is certainly 1. An LLR of  inf means the bit is 0.
* Output (of the decoder) is an int32_t array where each index is the error-corrected bit of the codeword at that index. 
* The termination condition of the iterative believe propagation algorithm is that the decoded codeword is a codeword of the parity check matrix.

