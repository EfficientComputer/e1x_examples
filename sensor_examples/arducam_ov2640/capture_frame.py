#!/usr/bin/env python3
"""Capture a JPEG frame streamed by the arducam single_frame example.

Reads the JPEG_BEGIN/<hex>/JPEG_END serial stream and writes a .jpg.
Requires pyserial (pip install pyserial). Flash/reset the board after starting.

Usage: python capture_frame.py -p /dev/ttyACM2 -o frame.jpg
"""
import argparse
import sys
import time

import serial  # pyserial


def main():
    ap = argparse.ArgumentParser(description="Capture an arducam JPEG frame over serial.")
    ap.add_argument("-p", "--port", default="/dev/ttyACM2", help="serial device")
    ap.add_argument("-b", "--baud", type=int, default=115200, help="baud rate")
    ap.add_argument("-o", "--out", default="frame.jpg", help="output JPEG path")
    ap.add_argument("-t", "--timeout", type=float, default=30.0,
                    help="seconds to wait for a frame")
    args = ap.parse_args()

    ser = serial.Serial(args.port, args.baud, timeout=1)
    print(f"listening on {args.port} @ {args.baud} — flash/reset the board now...")

    hex_chunks = []
    expected = None
    capturing = False
    deadline = time.time() + args.timeout

    while time.time() < deadline:
        line = ser.readline().decode("ascii", errors="replace").strip()
        if not line:
            continue

        if line.startswith("JPEG_BEGIN"):
            parts = line.split()
            expected = int(parts[1]) if len(parts) > 1 else None
            hex_chunks = []
            capturing = True
            print(f"frame start (expecting {expected} bytes)...")
            continue

        if line.startswith("JPEG_END"):
            if not capturing:
                continue
            data = bytes.fromhex("".join(hex_chunks))
            with open(args.out, "wb") as f:
                f.write(data)
            note = ""
            if expected is not None and len(data) != expected:
                note += f"  (WARNING: expected {expected})"
            note += "  [valid JPEG SOI]" if data[:2] == b"\xff\xd8" else "  [WARNING: no JPEG SOI]"
            print(f"wrote {len(data)} bytes to {args.out}{note}")
            return 0

        if capturing:
            hex_chunks.append("".join(c for c in line if c in "0123456789abcdefABCDEF"))
        else:
            print(line)  # echo other firmware output (banner, FIFO length, etc.)

    print("timed out waiting for a frame", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
