# Arducam OV2640 — E1x EVK

OV2640 + ArduChip (SPI FIFO) camera bring-up examples for the E1x EVK.

## Prerequisites

### I2C driver fix

The SDK I2C driver optimizes away buffer reads. Before building, add a compiler
barrier to `~/effcc/sdk/drivers/i2c/i2c.c` in `eff_i2c_read_raw()`:

```c
i2c_wait_int(i2c);

// Prevent compiler from optimizing away buffer reads
__asm__ volatile("" ::: "memory");

return 0;
```

### External level shifters

There is a communication issue with the TXB0108 on-board level shifters when interfacing with the ArduChip's SPI interface. We’ve seen the issue go away when using different external level shifters for the SPI interface. We specifically use the [4-channel I2C-safe Bi-directional Logic Level Converter BSS138](https://www.adafruit.com/product/757) level shifters from Adafruit.

## Wiring

The examples use EVK **bank 5** for SPI (ArduChip) and **bank 4 / I2C1** for I2C
(OV2640 sensor, SCCB address `0x30`).

| Signal | E1x EVK pin | Level shifting | Arducam pin |
|---|---|---|---|
| SPI CLK | DIGIO050 (bank 5, pin 0) | external BSS138 | SCK |
| SPI CS | DIGIO051 (bank 5, pin 1) | external BSS138 | CS |
| SPI MOSI | DIGIO052 (bank 5, pin 2) | external BSS138 | MOSI |
| SPI MISO | DIGIO053 (bank 5, pin 3) | external BSS138 | MISO |
| I2C SCL | DIGIO042 (bank 4, pin 2, I2C1) | external BSS138 | SCL |
| I2C SDA | DIGIO043 (bank 4, pin 3, I2C1) | external BSS138 | SDA |
| 3.3 V | 3V3 | — | VCC |
| GND | GND | — | GND |

```
E1x EVK (1.8 V IO)            BSS138 level shifter 1        Arducam Mini 2MP Plus
                              LV = 1.8 V      HV = 3.3 V
DIGIO050 (SPI CLK)  ────────▶ LV1 ─────────── HV1 ────────▶ SCK
DIGIO051 (SPI CS)   ────────▶ LV2 ─────────── HV2 ────────▶ CS
DIGIO052 (SPI MOSI) ────────▶ LV3 ─────────── HV3 ────────▶ MOSI
DIGIO053 (SPI MISO) ◀──────── LV4 ─────────── HV4 ◀──────── MISO

                              BSS138 level shifter 2
                              LV = 1.8 V      HV = 3.3 V
DIGIO042 (I2C1 SCL) ◀───────▶ LV1 ─────────── HV1 ◀───────▶ SCL
DIGIO043 (I2C1 SDA) ◀───────▶ LV2 ─────────── HV2 ◀───────▶ SDA

3.3 V ──▶ VCC (camera) and BSS138 HV;  1.8 V ──▶ BSS138 LV;  common GND
```

Notes:

- The four SPI lines need the external BSS138 shifters (see
  [External level shifters](#external-level-shifters) above); I2C works through
  the EVK's standard on-board level shifting, but we also level shift it via the
  external BSS138 level shifters in this example.
- SPI is mode 0, clocked at ~4.17 MHz (`clk_div = 5`), under the ArduChip's
  8 MHz max. I2C runs at 100 kHz.
- The E1x's digital IO is 1.8 V logic — never drive the EVK pins directly with
  the camera's 3.3 V signals.

## Bring-up order

Run the programs in this order — each adds one layer, so failures bisect cleanly:

| # | Program | Purpose |
|---|---------|---------|
| 1 | `arducam_id` | Bring up I2C; read the OV2640 product ID ("is it alive") |
| 2 | `spi_loopback` | Jumper DIGIO052 (MOSI) ↔ DIGIO053 (MISO) to verify SPI / level shifters |
| 3 | `cam_probe` | Bring up both: I2C PID on the sensor + SPI on the ArduChip |
| 4 | `fifo_bringup` | Capture + FIFO read path with no sensor config (a non-zero length proves the path) |
| 5 | `single_frame_jpeg_verify` | Full sensor config; read back 32 header bytes and verify the JPEG SOI |
| 6 | `single_frame` | Read a full frame from the camera and stream it over UART |

## Capturing an image (single_frame)

Start the host capture script *before* flashing/resetting the board:

```sh
python capture_frame.py -p /dev/ttyACM2 -o frame.jpg
```

Requires pyserial (`pip install pyserial`).
