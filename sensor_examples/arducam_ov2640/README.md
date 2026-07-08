# Arducam OV2640 — E1x EVK

OV2640 + ArduChip (SPI FIFO) camera bring-up examples for the E1x EVK.

## Prerequisite

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

There is a communication issue with the TXB0108 on-board level shifters when interfacing with the ArduChip's SPI interface. We’ve seen the issue go away when using different external level shifters for the SPI interface. We specifically use the 4-channel I2C-safe Bi-directional Logic Level Converter BSS138 level shifters from Adafruit.

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
