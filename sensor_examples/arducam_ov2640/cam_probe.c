//===-- cam_probe.c — Arducam OV2640+ArduChip probe: SPI test + I2C PID, E1x EVK ===//

#include <eff.h>
#include <eff/drivers/spi.h>
#include <eff/time.h>

#include <stdio.h>

// SPI / ArduChip: bank 5 (DIGIO050 CLK, 051 CS, 052 MOSI, 053 MISO).
#define CAM_SPI         SPI_5
#define CAM_SPI_PINMUX  PINMUX_5
#define CAM_SPI_CLKDIV  5  // 50MHz/(2*(5+1)) ~= 4.17MHz, under ArduCAM's 8MHz max

// I2C / OV2640: bank 4 pins 2/3 = I2C1 (DIGIO042 SCL, 043 SDA).
#define CAM_I2C         I2C_4_1
#define CAM_I2C_PINMUX  PINMUX_4
#define CAM_I2C_MODE    PINMUX_I2C0_I2C1
#define OV2640_ADDR     0x30

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

int main() {
    printf("\r\narducam / cam probe (bank5 SPI hw-CS, bank4 I2C, ext shifters)\r\n\r\n");

    // SPI ArduChip test register (ArduChip is SPI mode 0).
    eff_pinmux_set(CAM_SPI_PINMUX, PINMUX_SPI);
    eff_spi_cfg_t scfg = EFF_SPI_DEFAULTS;
    scfg.mode = SPI_MODE_0;
    scfg.xfer_mode = SPI_XFER_BIDIRECTIONAL;
    scfg.bus_size = SPI_BUS_SINGLE;
    scfg.clk_div = CAM_SPI_CLKDIV;
    eff_spi_init(CAM_SPI, &scfg);

    uint8_t tvals[3] = {0x55, 0xAA, 0x37};
    for (int i = 0; i < 3; i++) {
        arduchip_write(0x00, tvals[i]);
        uint8_t r = arduchip_read(0x00);
        printf("SPI test 0x%02x -> 0x%02x  %s\r\n", tvals[i], r,
               r == tvals[i] ? "OK" : "FAIL");
    }

    // I2C OV2640 product ID.
    eff_pinmux_set(CAM_I2C_PINMUX, CAM_I2C_MODE);
    if (eff_i2c_init(CAM_I2C, I2C_SPEED_100K)) {
        printf("I2C init failed\r\n");
        return -1;
    }

    uint8_t bb;
    bb = 0x01; eff_i2c_write(CAM_I2C, OV2640_ADDR, 0xFF, &bb, 1);  // sensor bank
    bb = 0x80; eff_i2c_write(CAM_I2C, OV2640_ADDR, 0x12, &bb, 1);  // soft reset
    sleep_ms(5);
    bb = 0x01; eff_i2c_write(CAM_I2C, OV2640_ADDR, 0xFF, &bb, 1);  // re-select bank
    sleep_us(500);

    // Settle between transactions: the OV2640 SCCB is slow to recover.
    uint8_t pidh = 0, pidl = 0;
    eff_i2c_read(CAM_I2C, OV2640_ADDR, 0x0A, &pidh, 1); sleep_us(500);
    eff_i2c_read(CAM_I2C, OV2640_ADDR, 0x0B, &pidl, 1);
    printf("PIDH = 0x%02x  (expect 0x26)\r\n", pidh);
    printf("PIDL = 0x%02x  (expect 0x41)\r\n", pidl);

    return 0;
}
