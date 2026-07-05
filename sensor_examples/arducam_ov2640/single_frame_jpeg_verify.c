//===-- single_frame_jpeg_verify.c — OV2640 VGA JPEG capture + header check, E1x EVK ===//
//
// Captures one 640x480 JPEG frame and verifies the SOI marker (0xFF 0xD8),
// dumping the first 32 bytes. Full-frame readout: single_frame.c.
//
//===----------------------------------------------------------------------===//

#include <eff.h>
#include <eff/drivers/spi.h>
#include <eff/time.h>

#include <stdio.h>

#include "ov2640_regs.h"

// SPI / ArduChip: bank 5 (DIGIO050 CLK, 051 CS, 052 MOSI, 053 MISO).
#define CAM_SPI         SPI_5
#define CAM_SPI_PINMUX  PINMUX_5
#define CAM_SPI_CLKDIV  5  // 50MHz/(2*(5+1)) ~= 4.17MHz, under ArduCAM's 8MHz max

// I2C / OV2640: bank 4 pins 2/3 = I2C1 (DIGIO042 SCL, 043 SDA).
#define CAM_I2C         I2C_4_1
#define CAM_I2C_PINMUX  PINMUX_4
#define CAM_I2C_MODE    PINMUX_I2C0_I2C1
#define OV2640_ADDR     0x30

// ArduChip registers.
#define ARDUCHIP_FIFO     0x04
#define FIFO_CLEAR_MASK   0x01
#define FIFO_START_MASK   0x02
#define ARDUCHIP_TRIG     0x41
#define CAP_DONE_MASK     0x08
#define FIFO_SIZE1        0x42
#define FIFO_SIZE2        0x43
#define FIFO_SIZE3        0x44
#define SINGLE_FIFO_READ  0x3D

#define READOUT_BYTES     32

// ArduChip access: write = reg|0x80, read = reg&0x7F (value lands on rx[1]).
static void arduchip_write(uint8_t reg, uint8_t val) {
    uint8_t tx[2] = {(uint8_t)(reg | 0x80), val};
    volatile uint8_t rx[2] = {0, 0};
    eff_spi_xfer(CAM_SPI, 0, 0, tx, 2, (uint8_t *)rx, 2);
}

static uint8_t arduchip_read(uint8_t reg) {
    uint8_t tx[2] = {(uint8_t)(reg & 0x7F), 0x00};
    volatile uint8_t rx[2] = {0, 0};
    eff_spi_xfer(CAM_SPI, 0, 0, tx, 2, (uint8_t *)rx, 2);
    return rx[1];
}

// Write an OV2640 register table (terminated by {0xff,0xff}).
static void ov2640_write_regs(const struct ov2640_reg *regs) {
    for (int i = 0; !(regs[i].reg == 0xff && regs[i].val == 0xff); i++) {
        uint8_t v = regs[i].val;
        eff_i2c_write(CAM_I2C, OV2640_ADDR, regs[i].reg, &v, 1);
        sleep_us(200);  // SCCB is slow; space writes out
    }
}

static uint32_t fifo_length(void) {
    uint32_t len = ((uint32_t)arduchip_read(FIFO_SIZE3) << 16) |
                   ((uint32_t)arduchip_read(FIFO_SIZE2) << 8) |
                   arduchip_read(FIFO_SIZE1);
    return len & 0x07ffff;
}

int main() {
    printf("\r\narducam / single frame JPEG verify (bank5 SPI, bank4 I2C)\r\n\r\n");

    // SPI to the ArduChip (mode 0).
    eff_pinmux_set(CAM_SPI_PINMUX, PINMUX_SPI);
    eff_spi_cfg_t scfg = EFF_SPI_DEFAULTS;
    scfg.mode = SPI_MODE_0;
    scfg.xfer_mode = SPI_XFER_BIDIRECTIONAL;
    scfg.bus_size = SPI_BUS_SINGLE;
    scfg.clk_div = CAM_SPI_CLKDIV;
    eff_spi_init(CAM_SPI, &scfg);

    // I2C to the OV2640.
    eff_pinmux_set(CAM_I2C_PINMUX, CAM_I2C_MODE);
    if (eff_i2c_init(CAM_I2C, I2C_SPEED_100K)) {
        printf("I2C init failed\r\n");
        return -1;
    }

    // Soft reset, then load the 640x480 (VGA) JPEG configuration.
    uint8_t bb;
    bb = 0x01; eff_i2c_write(CAM_I2C, OV2640_ADDR, 0xFF, &bb, 1);
    bb = 0x80; eff_i2c_write(CAM_I2C, OV2640_ADDR, 0x12, &bb, 1);
    sleep_ms(5);
    ov2640_write_regs(OV2640_JPEG_INIT);
    ov2640_write_regs(OV2640_YUV422);
    ov2640_write_regs(OV2640_JPEG);
    bb = 0x01; eff_i2c_write(CAM_I2C, OV2640_ADDR, 0xFF, &bb, 1);
    bb = 0x00; eff_i2c_write(CAM_I2C, OV2640_ADDR, 0x15, &bb, 1);
    ov2640_write_regs(OV2640_640x480_JPEG);
    sleep_ms(1000);  // let AEC/AGC/AWB converge to room light (early frames are dark)

    // Capture one frame: flush FIFO, start, wait for capture-done.
    arduchip_write(ARDUCHIP_FIFO, FIFO_CLEAR_MASK);
    arduchip_write(ARDUCHIP_FIFO, FIFO_START_MASK);

    int timeout = 1000000;
    while (!(arduchip_read(ARDUCHIP_TRIG) & CAP_DONE_MASK) && --timeout) {}
    if (timeout == 0) {
        printf("capture timed out (no CAP_DONE)\r\n");
        return -1;
    }

    uint32_t len = fifo_length();
    printf("capture done, FIFO length = %lu bytes\r\n", (unsigned long)len);
    if (len == 0) {
        printf("FIFO empty\r\n");
        return 0;
    }

    // Read out via single FIFO read (no dummy discard — first byte is the SOI).
    uint8_t b0 = arduchip_read(SINGLE_FIFO_READ);
    uint8_t b1 = arduchip_read(SINGLE_FIFO_READ);
    printf("SOI = 0x%02x%02x  %s\r\n", b0, b1,
           (b0 == 0xff && b1 == 0xd8) ? "valid JPEG" : "NOT a JPEG SOI");

    printf("first %d bytes: %02x %02x", READOUT_BYTES, b0, b1);
    for (int i = 2; i < READOUT_BYTES; i++) {
        printf(" %02x", arduchip_read(SINGLE_FIFO_READ));
    }
    printf("\r\n");

    return 0;
}
