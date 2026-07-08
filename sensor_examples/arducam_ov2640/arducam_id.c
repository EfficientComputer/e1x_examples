//===-- arducam_id.c — Arducam OV2640 I2C "is it alive" PID check, E1x EVK ===//

#include <eff.h>

#include <stdio.h>

// I2C / OV2640: bank 4 pins 2/3 = I2C1 (DIGIO042 SCL, 043 SDA), SW10 ON.
#define CAM_I2C         I2C_4_1
#define CAM_PINMUX      PINMUX_4
#define CAM_PINMUX_MODE PINMUX_I2C0_I2C1
#define OV2640_ADDR     0x30

#define OV2640_REG_BANK_SEL 0xFF
#define OV2640_BANK_SENSOR  0x01
#define OV2640_REG_PIDH     0x0A  // expect 0x26
#define OV2640_REG_PIDL     0x0B  // expect 0x41

int main() {
    printf("\r\narducam / OV2640 id check\r\n\r\n");

    eff_pinmux_set(CAM_PINMUX, CAM_PINMUX_MODE);
    eff_i2c_init(CAM_I2C, I2C_SPEED_100K);

    uint8_t bank = OV2640_BANK_SENSOR;  // PIDH/PIDL live in the sensor bank
    eff_i2c_write(CAM_I2C, OV2640_ADDR, OV2640_REG_BANK_SEL, &bank, 1);

    uint8_t pidh = 0, pidl = 0;
    eff_i2c_read(CAM_I2C, OV2640_ADDR, OV2640_REG_PIDH, &pidh, 1);
    eff_i2c_read(CAM_I2C, OV2640_ADDR, OV2640_REG_PIDL, &pidl, 1);
    printf("PIDH = 0x%02x  (expect 0x26)\r\n", pidh);
    printf("PIDL = 0x%02x  (expect 0x41)\r\n", pidl);

    return 0;
}
