#!/usr/bin/env python
"""
gen_mel_filterbank.py - Generate a C header for the mel filterbank.

Outputs a MelFilter struct definition, a flat int16 weights array, and a
static array of MelFilter structs (one per mel bin).

Usage
-----
  python gen_mel_filterbank.py                        # defaults, to stdout
  python gen_mel_filterbank.py -o mel_filterbank.h
  python gen_mel_filterbank.py --num_mel_bins 40 --window_size_samples 1024 --sample_rate 16000
"""

import argparse
import subprocess
import sys

import numpy as np
import tensorflow as tf


def gen_mel_filterbank_c(
    num_mel_bins=40, num_spectrogram_bins=513, sample_rate=16000, out=None
):
    """Print C source for the mel filterbank to *out* (default: stdout)."""
    if out is None:
        out = sys.stdout

    lower_edge_hertz = 0.0
    upper_edge_hertz = float(sample_rate) / 2.0

    weight_matrix = tf.signal.linear_to_mel_weight_matrix(
        num_mel_bins=num_mel_bins,
        num_spectrogram_bins=num_spectrogram_bins,
        sample_rate=sample_rate,
        lower_edge_hertz=lower_edge_hertz,
        upper_edge_hertz=upper_edge_hertz,
    )

    mel_matrix_int16 = np.round(weight_matrix.numpy() * np.iinfo(np.int16).max).astype(
        np.int16
    )

    # Build per-filter (start_bin, weights) from each column of the matrix.
    filters = []
    all_weights = []
    weight_offsets = []
    for f in range(num_mel_bins):
        col = mel_matrix_int16[:, f]
        nonzero_rows = np.flatnonzero(col)
        if len(nonzero_rows) == 0:
            start_bin, weights = 0, []
        else:
            start_bin = int(nonzero_rows[0])
            end_bin = int(nonzero_rows[-1])
            weights = col[start_bin : end_bin + 1].tolist()
        weight_offsets.append(len(all_weights))
        all_weights.extend(weights)
        filters.append((start_bin, weights))

    total_weights = len(all_weights)

    def w(s=""):
        out.write(s + "\n")

    w("#pragma once")
    w("#include <stdint.h>")
    w()
    w(
        f"/* Mel filterbank: {num_mel_bins} triangular filters over "
        f"{num_spectrogram_bins} spectrogram bins,"
    )
    w(f"   sample_rate={sample_rate} Hz, weights normalized to int16. */")
    w(
        f"/* Total stored weights: {total_weights} "
        f"(vs {num_spectrogram_bins * num_mel_bins} dense) */"
    )
    w()
    w(f"#define NUM_MEL_FILTERS          {num_mel_bins}")
    w(f"#define MEL_FILTER_TOTAL_WEIGHTS {total_weights}")
    w()

    # Flat weights array
    w(f"static const int16_t mel_filter_weights[{total_weights}] = {{")
    chunk = 16
    for i in range(0, total_weights, chunk):
        seg = all_weights[i : i + chunk]
        comma = "," if i + chunk < total_weights else ""
        w(f'    {", ".join(str(x) for x in seg)}{comma}')
    w("};")
    w()

    # Struct definition + array
    w("typedef struct {")
    w("    int16_t location;")
    w("    int16_t length;")
    w("    int16_t offset; /* index into mel_filter_weights[] */")
    w("} MelFilter;")
    w()
    w(f"static const MelFilter mel_filters[{num_mel_bins}] = {{")
    for f, (start_bin, weights) in enumerate(filters):
        comma = "," if f < num_mel_bins - 1 else ""
        w(
            f"    {{ {start_bin:3d}, {len(weights):3d}, {weight_offsets[f]:4d} }}{comma}"
            f"  /* filter {f} */"
        )
    w("};")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--num_mel_bins", type=int, default=40)
    parser.add_argument(
        "--window_size_samples",
        type=int,
        default=1024,
        help="STFT window length in samples (default 1024). "
        "num_spectrogram_bins = window_size_samples // 2 + 1",
    )
    parser.add_argument("--sample_rate", type=int, default=16000)
    parser.add_argument(
        "-o", "--output", default=None, help="Output file path (default: stdout)"
    )
    args = parser.parse_args()

    num_spectrogram_bins = args.window_size_samples // 2 + 1

    if args.output:
        with open(args.output, "w") as f:
            gen_mel_filterbank_c(
                args.num_mel_bins, num_spectrogram_bins, args.sample_rate, out=f
            )
        subprocess.run(["clang-format", "-style=file", "-i", args.output], check=True)
        print(f"Wrote {args.output}")
    else:
        gen_mel_filterbank_c(args.num_mel_bins, num_spectrogram_bins, args.sample_rate)


if __name__ == "__main__":
    main()
