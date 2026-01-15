/*
 * ADXL372 Accelerometer Example for E1x
 * 
 * Minimal SPI interface to read accelerometer data.
 * No floating point operations used.
 */

#include <eff.h>
#include <eff/drivers/spi.h>
#include <stdio.h>
#include <stdint.h>

/* E1x Configuration */
#define SPI_PORT                SPI_2
#define PINMUX_PORT             PINMUX_2

/* ADXL372 Register Addresses */
#define ADXL372_DEVID           0x00
#define ADXL372_DEVID_MST       0x01
#define ADXL372_PARTID          0x02
#define ADXL372_STATUS_1        0x04
#define ADXL372_X_DATA_H        0x08
#define ADXL372_X_DATA_L        0x09
#define ADXL372_Y_DATA_H        0x0A
#define ADXL372_Y_DATA_L        0x0B
#define ADXL372_Z_DATA_H        0x0C
#define ADXL372_Z_DATA_L        0x0D
#define ADXL372_TIMING          0x3D
#define ADXL372_MEASURE         0x3E
#define ADXL372_POWER_CTL       0x3F
#define ADXL372_RESET           0x41

/* Expected ID values */
#define ADXL372_DEVID_VAL       0xAD
#define ADXL372_PARTID_VAL      0xFA

/* SPI read/write bits */
#define ADXL372_READ            0x01
#define ADXL372_WRITE           0x00

/* Operating modes */
#define ADXL372_STANDBY         0x00
#define ADXL372_FULL_BW         0x03

/* Convert 12-bit two's complement to signed 16-bit */
static int16_t twos_complement_12bit(uint16_t val) {
    if (val & 0x800) {
        return (int16_t)(val | 0xF000);
    }
    return (int16_t)val;
}

/* Read a single register */
static uint8_t adxl372_read_reg(uint8_t reg) {
    uint8_t tx = (reg << 1) | ADXL372_READ;
    uint8_t rx = 0;
    eff_spi_xfer(SPI_PORT, 0, 0, &tx, 1, &rx, 1);
    return rx;
}

/* Write a single register */
static void adxl372_write_reg(uint8_t reg, uint8_t val) {
    uint8_t tx[2];
    tx[0] = (reg << 1) | ADXL372_WRITE;
    tx[1] = val;
    eff_spi_xfer(SPI_PORT, 0, 0, tx, 2, NULL, 0);
}

/* Read multiple registers */
static void adxl372_read_multi(uint8_t reg, uint8_t *buf, uint8_t len) {
    uint8_t tx = (reg << 1) | ADXL372_READ;
    eff_spi_xfer(SPI_PORT, 0, 0, &tx, 1, buf, len);
}

/* Initialize ADXL372 */
static int adxl372_init(void) {
    uint8_t dev_id, part_id;
    
    /* Reset the device */
    adxl372_write_reg(ADXL372_RESET, 0x52);
    sleep_ms(100);  /* Wait for reset to complete */
    
    /* First read after reset is unreliable - discard it */
    (void)adxl372_read_reg(ADXL372_DEVID);
    
    /* Verify device ID */
    dev_id = adxl372_read_reg(ADXL372_DEVID);
    part_id = adxl372_read_reg(ADXL372_PARTID);
    
    if (dev_id != ADXL372_DEVID_VAL || part_id != ADXL372_PARTID_VAL) {
        printf("Error: Invalid device ID (0x%02X, 0x%02X)\r\n", dev_id, part_id);
        return -1;
    }
    
    printf("ADXL372 detected (ID: 0x%02X, Part: 0x%02X)\r\n", dev_id, part_id);
    
    /* Configure: ODR 400Hz, bandwidth 200Hz */
    adxl372_write_reg(ADXL372_TIMING, 0x00);
    adxl372_write_reg(ADXL372_MEASURE, 0x00);
    
    /* Set to full bandwidth measurement mode */
    adxl372_write_reg(ADXL372_POWER_CTL, ADXL372_FULL_BW);
    
    sleep_ms(10);
    return 0;
}

/* Read acceleration data (returns raw 12-bit values) */
static void adxl372_read_accel(int16_t *x, int16_t *y, int16_t *z) {
    uint8_t buf[6];
    uint16_t raw_x, raw_y, raw_z;
    
    /* Wait for data ready */
    while (!(adxl372_read_reg(ADXL372_STATUS_1) & 0x01)) {
        /* spin */
    }
    
    /* Read all 6 bytes at once */
    adxl372_read_multi(ADXL372_X_DATA_H, buf, 6);
    
    /* Combine bytes: data is 12-bit, upper 8 bits in H register, lower 4 in L */
    raw_x = ((uint16_t)buf[0] << 4) | (buf[1] >> 4);
    raw_y = ((uint16_t)buf[2] << 4) | (buf[3] >> 4);
    raw_z = ((uint16_t)buf[4] << 4) | (buf[5] >> 4);
    
    /* Convert to signed */
    *x = twos_complement_12bit(raw_x);
    *y = twos_complement_12bit(raw_y);
    *z = twos_complement_12bit(raw_z);
}

int main(void) {
    int16_t x, y, z;
    eff_spi_cfg_t spi_cfg = EFF_SPI_DEFAULTS;
    
    printf("ADXL372 Accelerometer Example\r\n");
    
    /* Configure pinmux for SPI */
    eff_pinmux_set(PINMUX_PORT, PINMUX_SPI);
    
    /* Initialize SPI: write then read mode, single bus */
    spi_cfg.xfer_mode = SPI_XFER_WRITE_READ;
    spi_cfg.bus_size = SPI_BUS_SINGLE;
    spi_cfg.clk_div = 4;  /* Adjust for desired SPI clock speed */
    
    if (eff_spi_init(SPI_PORT, &spi_cfg) != 0) {
        printf("Error: SPI init failed\r\n");
        return -1;
    }
    
    /* Initialize ADXL372 */
    if (adxl372_init() != 0) {
        return -1;
    }
    
    printf("Reading accelerometer data...\r\n");
    printf("(Raw 12-bit values, scale: 100 mg/LSB)\r\n\r\n");
    
    /* Main loop: read and print acceleration data */
    while (1) {
        adxl372_read_accel(&x, &y, &z);
        
        /* Print raw values (multiply by 100 to get milli-g) */
        printf("X:%3d\t  Y:%3d\t  Z:%3d\r\n", x, y, z);
        
        sleep_ms(1000);
    }
    
    return 0;
}

