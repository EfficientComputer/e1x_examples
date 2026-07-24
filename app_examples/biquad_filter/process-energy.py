import argparse
import csv
import sys

# The EVK firmware prints this header once at startup. If the capture began
# after it scrolled past (a common minicom race), fall back to this so the
# capture is still usable without re-capturing. Keep in sync with the EVK data
# format; see docs.efficient.computer/evaluation-kit#collecting-energy-data.
DEFAULT_EVK_HEADER = (
    "timestamp(us),current_sys(mA),voltage_sys(mV),power_sys(mW),"
    "current_1v8(mA),voltage_1v8(mV),power_1v8(mW),"
    "current_e1x(mA),voltage_e1x(mV),power_e1x(mW),"
    "current_var(mA),voltage_var(mV),power_var(mW),aon4,aon5"
).split(",")


def read_capture(path):
    header = None
    rows = []

    with open(path, "r", encoding="utf-8", errors="ignore") as capture:
        for line in capture:
            line = line.strip()
            if not line:
                continue
            if line.startswith("timestamp(us)"):
                header = line.split(",")
                continue
            if line[0].isdigit():
                rows.append(line)

    if header is None:
        # No header in the capture (e.g. it started mid-stream). Fall back to
        # the known EVK column layout, but only for rows whose field count
        # matches so a truncated first/last line cannot misalign the columns.
        expected = len(DEFAULT_EVK_HEADER)
        rows = [row for row in rows if len(row.split(",")) == expected]
        if not rows:
            raise RuntimeError("No EVK CSV header found and no usable data rows in capture")
        header = DEFAULT_EVK_HEADER
        print(
            "warning: no header in capture; assuming default EVK column layout",
            file=sys.stderr,
        )

    return list(csv.DictReader(rows, fieldnames=header))


def signal_value(row, signal):
    return int(float(row[signal]))


def find_low_region(samples, signal):
    armed = False
    previous = None
    start = None

    for index, row in enumerate(samples):
        value = signal_value(row, signal)
        if value == 1:
            armed = True
        if armed and previous == 1 and value == 0:
            start = index
            break
        previous = value

    if start is None:
        raise RuntimeError(f"No falling edge found on {signal}")

    end = start
    while end + 1 < len(samples) and signal_value(samples[end + 1], signal) == 0:
        end += 1

    return samples[start : end + 1]


def mean(rows, column):
    return sum(float(row[column]) for row in rows) / len(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("capture")
    parser.add_argument("--iters", type=int, required=True)
    parser.add_argument("--signal", default="aon4")
    parser.add_argument("--energy-column", default="power_e1x(mW)")
    args = parser.parse_args()

    samples = read_capture(args.capture)
    run = find_low_region(samples, args.signal)

    start_us = float(run[0]["timestamp(us)"])
    end_us = float(run[-1]["timestamp(us)"])
    runtime_ms = (end_us - start_us) / 1000.0

    mean_power_mw = mean(run, args.energy_column)
    total_energy_uj = mean_power_mw * runtime_ms

    print(f"Runtime: {runtime_ms:.3f} ms")
    print(f"Iterations: {args.iters}")
    print(f"Time per iteration: {runtime_ms / args.iters:.6f} ms")
    print(f"Mean {args.energy_column}: {mean_power_mw:.6f} mW")
    print(f"Energy per iteration: {total_energy_uj / args.iters:.6f} uJ")

    for column in ("power_sys(mW)", "power_1v8(mW)", "power_e1x(mW)", "power_var(mW)"):
        if column in run[0]:
            print(f"Mean {column}: {mean(run, column):.6f} mW")


if __name__ == "__main__":
    main()
