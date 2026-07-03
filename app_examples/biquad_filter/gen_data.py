import argparse
import math

def write_array(f, arr, per_line=8):
    for i, v in enumerate(arr):
        if i != 0:
            f.write(", ")

        if i % per_line == per_line - 1:
            f.write("\n\t")

        f.write(f"{v}")


def gen_sample_input(size, f1, f2, scale=15000, phase_shift=0):
    ret = []

    s1_phase_shift = phase_shift
    s2_phase_shift = phase_shift
    if f1 == 0:
        s1_phase_shift = 0
    if f2 == 0:
        s2_phase_shift = 0

    for i in range(size):
        s1 = math.sin(i / size * f1 * 2 * math.pi + s1_phase_shift) * scale
        s2 = math.sin(i / size * f2 * 2 * math.pi + s2_phase_shift) * scale
        ret.append(int(s1 + s2))

    return ret


def generate_sample(path, size, f1, f2, low_gain, phase_shift):
    with open(path, "w") as f:
        f.write("#include <stdint.h>\n\n")
        f.write(f"const int16_t sample_input[{size}] = {{")
        write_array(f, gen_sample_input(size, f1, f2))
        f.write("};\n\n")
        f.write(f"const int16_t expected_output[{size}] = {{")
        low = min(f1, f2)
        write_array(f, gen_sample_input(size, low, 0, scale=15000 * low_gain, phase_shift=phase_shift))
        f.write("};\n\n")

    # Make a plot to visualize the input
    try:
        import matplotlib.pyplot as plt

        plt.plot(gen_sample_input(size, f1, f2)[:512])
        plt.title("Sample Input")
        plt.xlabel("Sample")
        plt.ylabel("Amplitude")
        plt.grid()
        plt.savefig(path + ".png")
        plt.close()
    except ImportError:
        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--length", default=512, type=int, help="Number of samples"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Output file path"
    )
    parser.add_argument(
        "--freq1", required=True, type=int, help="Wave one frequency relative to samples"
    )
    parser.add_argument(
        "--freq2", required=True, type=int, help="Wave two frequency relative to samples"
    )
    parser.add_argument(
        "--low-gain", default=1, type=float, help="Low passed wave gain"
    )
    parser.add_argument(
        "--phase-shift", default=1, type=float, help="Low passed wave phase shift"
    )
    args = parser.parse_args()

    generate_sample(args.output, args.length, args.freq1, args.freq2, args.low_gain, args.phase_shift)
