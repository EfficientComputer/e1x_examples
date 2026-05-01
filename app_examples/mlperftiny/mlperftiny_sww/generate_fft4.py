import argparse
import math
import subprocess
import sys


def twiddle_size(n, radix):
    if n <= radix:
        return radix
    return n + twiddle_size(n / radix, radix)


def gen_twiddle_schedule(s, n, radix):
    if n <= radix:
        return [int(s)]

    rest = gen_twiddle_schedule(s + n, n / radix, radix)
    rest.extend([int(s)])
    return rest


def gen_twiddle(idx, n, inverse):
    phase = -2 * math.pi * idx / n
    if inverse:
        phase = -phase

    return (
        math.floor(0.5 + 32767 * math.cos(phase)),
        math.floor(0.5 + 32767 * math.sin(phase)),
    )


def gen_twiddles(n, inverse, radix):
    twiddles = []

    for i in range(n):
        twiddle = gen_twiddle(i, n, inverse)
        twiddles.append(twiddle[0])
        twiddles.append(twiddle[1])

    j = int(n / radix)
    while j >= radix:
        for i in range(j):
            twiddle = gen_twiddle(i * n / j, n, inverse)
            twiddles.append(twiddle[0])
            twiddles.append(twiddle[1])

        j = int(j / radix)

    return twiddles


def write_array(f, arr, per_line=8):
    for i, v in enumerate(arr):
        if i != 0:
            f.write(", ")

        if i % per_line == per_line - 1:
            f.write("\n\t")

        f.write(f"{v}")


def generate_fft4(size, path, inverse, radix):
    with open(path, "w") as f:
        f.write("#include <stdint.h>\n\n")
        f.write(f"#define FFT_SIZE {size}\n")
        f.write(f"#define FFT_IS_INVERSE {1 if inverse else 0}\n\n")
        f.write(f"#ifdef DEFINE_TWIDDLES\n")
        f.write(f"const int32_t effFFTDummyWord = 0xEFF;\n")
        f.write(
            f"const int16_t __attribute__((aligned(4))) twiddles[{int(twiddle_size(size, radix)) * 2}] = {{\n\t"
        )

        write_array(f, gen_twiddles(size, inverse, radix))

        f.write("};\n\n")

        tw_schedule = gen_twiddle_schedule(0, size, radix)
        f.write(f"const int twiddleSchedule[{len(tw_schedule)}] = {{\n\t")
        write_array(f, tw_schedule)

        f.write("};\n\n")
        f.write("#endif\n\n")

    subprocess.run(["clang-format", "-style=file", "-i", path], check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--size", default=256, type=int, help="Size of fft. Must be power of 4"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Output file path"
    )
    parser.add_argument(
        "-x", "--inverse", default=False, type=bool, help="Whether the FFT is inverted"
    )
    parser.add_argument(
        "-c",
        "--chip",
        required=True,
        type=str,
        help="One of 'e0' or 'e1x' or 'e1' or 'e2'",
    )
    args = parser.parse_args()

    if args.chip in ["e1", "e1x", "e2"]:
        radix = 4
    elif args.chip == "e0":
        radix = 2
    else:
        print(f"unknown chip target [{args.chip}]", file=sys.stderr)

    generate_fft4(args.size, args.output, args.inverse, radix)
