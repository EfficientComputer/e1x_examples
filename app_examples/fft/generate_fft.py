import argparse
import math
import sys


def twiddle_size(n, radix):
    if n <= radix:
        return radix
    return n + twiddle_size(n / radix, radix)


def gen_twiddle_schedule(s, n, radix):
    if n <= radix:
        return [s]

    rest = gen_twiddle_schedule(s + n, n / radix, radix)
    rest.extend([s])
    return rest


def gen_twiddle(idx, n, inverse):
    phase = -2 * math.pi * idx / n
    if inverse:
        phase = -phase

    return (
        math.floor(0.5 + 32767 * math.cos(phase)),
        math.floor(0.5 + 32767 * math.sin(phase)),
    )


def gen_twiddles(n, inverse, radix, transposed=False):
    twiddles = []

    def emit(pairs):
        if transposed:
            for k in range(0, len(pairs), 2):
                r0, i0 = pairs[k]
                r1, i1 = pairs[k + 1]
                twiddles.extend([r0, r1, i0, i1])
        else:
            for r, i in pairs:
                twiddles.extend([r, i])

    emit([gen_twiddle(i, n, inverse) for i in range(n)])

    j = int(n / radix)
    while j >= radix:
        emit([gen_twiddle(i * n / j, n, inverse) for i in range(j)])
        j = int(j / radix)

    return twiddles


def gen_sample_input(size):
    ret = []

    for i in range(size):
        s100 = math.sin(i / size * 100) * 10000
        s30 = math.sin(i / size * 30) * 10000
        s1000 = math.sin(i / size * 1000) * 10000
        ret.append(int(s100 + s30 + s1000))

    return ret


def write_array(f, arr, per_line=8):
    for i, v in enumerate(arr):
        if i != 0:
            f.write(", ")

        if i % per_line == per_line - 1:
            f.write("\n\t")

        f.write(f"{v}")


def _int16(x):
    x = int(x) & 0xFFFF
    return x - 0x10000 if x >= 0x8000 else x


def _sround(x):
    return _int16((int(x) + (1 << 14)) >> 15)


def _c_fixdiv(r, i, divisor):
    scale = 32767 // divisor
    return _sround(r * scale), _sround(i * scale)


def _c_mul(ar, ai, br, bi):
    return _sround(ar * br - ai * bi), _sround(ar * bi + ai * br)


def _unpack_u32(v):
    v = int(v) & 0xFFFFFFFF
    r = v & 0xFFFF
    im = (v >> 16) & 0xFFFF
    return _int16(r), _int16(im)


def _digit_rev_radix4(i, log2size):
    rev32 = int(f"{i:032b}"[::-1], 2)
    a = rev32 & 0xAAAAAAAA
    b = rev32 & 0x55555555
    return (a >> (32 - log2size + 1)) | (b >> (32 - log2size - 1))


def _digit_rev_radix2(i, log2size):
    rev32 = int(f"{i:032b}"[::-1], 2)
    return rev32 >> (32 - log2size)


def simulate_fft(size, inverse, radix):
    import math

    log2size = int(math.log2(size) + 0.5)
    sample = gen_sample_input(size)
    tw_flat = gen_twiddles(size, inverse, radix)
    twiddles = [(tw_flat[2 * k], tw_flat[2 * k + 1]) for k in range(len(tw_flat) // 2)]
    sched = gen_twiddle_schedule(0, size, radix)

    data = [
        _unpack_u32(
            sample[
                (
                    _digit_rev_radix4(i, log2size)
                    if radix == 4
                    else _digit_rev_radix2(i, log2size)
                )
            ]
        )
        for i in range(size)
    ]

    layer, m = 0, 1
    while m * radix <= size:
        ts = int(sched[layer])
        for seg in range(size // (m * radix)):
            base = seg * m * radix
            for j in range(m):
                if radix == 4:
                    f0r, f0i = _c_fixdiv(*data[base + j], 4)
                    f1r, f1i = _c_fixdiv(*data[base + j + m], 4)
                    f2r, f2i = _c_fixdiv(*data[base + j + 2 * m], 4)
                    f3r, f3i = _c_fixdiv(*data[base + j + 3 * m], 4)
                    s0r, s0i = _c_mul(f1r, f1i, *twiddles[ts + j])
                    s1r, s1i = _c_mul(f2r, f2i, *twiddles[ts + 2 * j])
                    s2r, s2i = _c_mul(f3r, f3i, *twiddles[ts + 3 * j])
                    s5r, s5i = _int16(f0r - s1r), _int16(f0i - s1i)
                    f0r, f0i = _int16(f0r + s1r), _int16(f0i + s1i)
                    s3r, s3i = _int16(s0r + s2r), _int16(s0i + s2i)
                    s4r, s4i = _int16(s0r - s2r), _int16(s0i - s2i)
                    data[base + j] = (_int16(f0r + s3r), _int16(f0i + s3i))
                    data[base + j + 2 * m] = (_int16(f0r - s3r), _int16(f0i - s3i))
                    if inverse:
                        data[base + j + m] = (_int16(s5r - s4i), _int16(s5i + s4r))
                        data[base + j + 3 * m] = (_int16(s5r + s4i), _int16(s5i - s4r))
                    else:
                        data[base + j + m] = (_int16(s5r + s4i), _int16(s5i - s4r))
                        data[base + j + 3 * m] = (_int16(s5r - s4i), _int16(s5i + s4r))
                else:
                    f0r, f0i = _c_fixdiv(*data[base + j], 2)
                    f1r, f1i = _c_fixdiv(*data[base + j + m], 2)
                    s0r, s0i = _c_mul(f1r, f1i, *twiddles[ts + j])
                    data[base + j + m] = (_int16(f0r - s0r), _int16(f0i - s0i))
                    data[base + j] = (_int16(f0r + s0r), _int16(f0i + s0i))
        m *= radix
        layer += 1

    return data


def validate(size, inverse, radix):
    import numpy as np

    result = simulate_fft(size, inverse, radix)
    fixed_R = np.array([r for r, _ in result], dtype=float)
    fixed_I = np.array([i for _, i in result], dtype=float)
    sample = gen_sample_input(size)
    cpx_in = np.array([_unpack_u32(v)[0] + 1j * _unpack_u32(v)[1] for v in sample])
    ref = np.fft.fft(cpx_in) / size
    err_R = fixed_R - ref.real
    err_I = fixed_I - ref.imag
    all_err = np.concatenate([err_R, err_I])
    max_err = np.max(np.abs(all_err))
    rms_err = np.sqrt(np.mean(all_err**2))
    sig_pow = np.mean(ref.real**2 + ref.imag**2)
    snr = 10 * np.log10(sig_pow / np.mean(all_err**2)) if rms_err > 0 else float("inf")
    return max_err, rms_err, snr


def generate_expected(size, path, inverse, radix):
    result = simulate_fft(size, inverse, radix)
    with open(path, "w") as f:
        f.write("#include <stdint.h>\n\n")
        f.write(f"int32_t expectedR[{size}] = {{\n\t")
        write_array(f, [r for r, _ in result])
        f.write("};\n\n")
        f.write(f"int32_t expectedI[{size}] = {{\n\t")
        write_array(f, [i for _, i in result])
        f.write("};\n")


def generate_fft(size, path, inverse, radix, transposed=False):
    with open(path, "w") as f:
        f.write("#include <stdint.h>\n\n")
        f.write(f"#define FFT_IS_INVERSE {1 if inverse else 0}\n\n")
        f.write(
            f"const int16_t _Alignas(4) twiddles[{int(twiddle_size(size, radix)) * 2}] = {{\n\t"
        )

        write_array(f, gen_twiddles(size, inverse, radix, transposed))

        f.write("};\n\n")

        tw_schedule = gen_twiddle_schedule(0, size, radix)
        f.write(f"const int twiddleSchedule[{len(tw_schedule)}] = {{\n\t")
        write_array(f, tw_schedule)

        f.write("};\n\n")

        f.write(f"const uint32_t sample_input[{size}] = {{")
        write_array(f, gen_sample_input(size))
        f.write("};\n\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--size", default=256, type=int, help="Size of fft. Must be power of 4"
    )
    parser.add_argument(
        "-o", "--output", default=None, type=str, help="Output file path"
    )
    parser.add_argument(
        "-x", "--inverse", default=False, type=bool, help="Whether the FFT is inverted"
    )
    parser.add_argument(
        "-t",
        "--transposed",
        default=False,
        action="store_true",
        help="Output twiddles with 2x2 transposed layout: [r0,r1,i0,i1] per pair instead of [r0,i0,r1,i1]",
    )
    parser.add_argument(
        "-e",
        "--expected",
        default=False,
        action="store_true",
        help="Generate expected output arrays (expectedR, expectedI) instead of twiddle tables",
    )
    parser.add_argument(
        "--validate",
        default=False,
        action="store_true",
        help="Compare fixed-point simulation against numpy FFT reference and print error metrics",
    )
    args = parser.parse_args()

    if args.validate:
        max_err, rms_err, snr = validate(args.size, args.inverse, radix)
        print(
            f"{'size':>6}  {'chip':<14}  {'max_err':>8}  {'rms_err':>8}  {'SNR(dB)':>8}"
        )
        print(
            f"{args.size:>6}  {'e1x':<14}  {max_err:>8.2f}  {rms_err:>8.3f}  {snr:>8.1f}"
        )
    elif args.expected:
        if args.output is None:
            print(
                "error: -o/--output is required when not using --validate",
                file=sys.stderr,
            )
            sys.exit(1)

        radix = 4
        generate_expected(args.size, args.output, args.inverse, radix)
    else:
        if args.output is None:
            print(
                "error: -o/--output is required when not using --validate",
                file=sys.stderr,
            )
            sys.exit(1)

        radix = 4
        generate_fft(args.size, args.output, args.inverse, radix, args.transposed)
