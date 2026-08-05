/**
 ******************************************************************************
 *
 * @file neo_m9n.h
 *
 * @brief u-blox NEO-M9N GNSS receiver driver — public API
 *
 * Interface: I2C (u-blox DDC)
 * I2C address: 0x42 (default, configurable via UBX-CFG-I2C)
 * Host: Efficient Computer E1x rev3 EVK
 *
 * The NEO-M9N streams NMEA sentences via a FIFO accessible through three
 * registers:
 *   0xFD  high byte of bytes available in output FIFO
 *   0xFE  low  byte of bytes available in output FIFO
 *   0xFF  data stream — repeated reads advance through the FIFO
 *
 * NOTE: unlike the Bosch BMI323/BMM350, the NEO-M9N does NOT insert dummy
 * bytes before its payload.  Reads from 0xFD/0xFF return actual data in
 * buf[0] onward.
 *
 ******************************************************************************
 */

#ifndef NEO_M9N_H
#define NEO_M9N_H

#include <eff.h>
#include <eff/drivers/i2c.h>
#include <stdint.h>

/*
 * MACROS
 ******************************************************************************
 */
#define NEO_M9N_I2C_ADDR    0x42u

/* Parsed subset of UBX-NAV-PVT (see neo_m9n_poll_pvt). */
typedef struct {
    uint8_t fix_type;    /* 0 = none, 2 = 2D, 3 = 3D */
    uint8_t num_sv;      /* satellites used in solution */
    int32_t lon_1e7;     /* longitude, 1e-7 deg */
    int32_t lat_1e7;     /* latitude,  1e-7 deg */
    int32_t height_mm;   /* height above ellipsoid, mm */
} neo_m9n_pvt_t;

/*
 * PUBLIC API
 ******************************************************************************
 */

/**
 * @brief  Initialise pinmux and I2C peripheral for the NEO-M9N.
 *         Call once at startup before any other neo_m9n_* function.
 * @return 0 on success, -1 on I2C init error.
 */
int8_t neo_m9n_init(void);

/**
 * @brief  Query how many bytes are waiting in the NEO-M9N output FIFO.
 * @return Number of bytes available (0 = FIFO empty), or -1 on I2C error
 *         (device not present / no ACK at 0x42).
 */
int32_t neo_m9n_bytes_available(void);

/**
 * @brief  Send a UBX-CFG-VALSET command to enable the antenna supply voltage
 *         on the RF pin (CFG-HW-ANT_CFG_VOLTCTRL = 1).  Required to power
 *         active (LNA) antennas.  Call once after neo_m9n_init().
 * @return 0 on success, -1 on I2C write error.
 */
int8_t neo_m9n_enable_active_antenna(void);

/**
 * @brief  Send a UBX-CFG-VALSET command to enable NMEA output on the I2C DDC
 *         port (CFG-I2COUTPROT-NMEA = 1).  Call once after neo_m9n_init().
 *         Required because the M9 platform does not enable I2C NMEA output by
 *         default (CFG-PRT is deprecated on M9 — use this instead).
 * @return 0 on success, -1 on I2C write error.
 */
int8_t neo_m9n_enable_nmea_i2c(void);

/**
 * @brief  Disable NMEA output on the I2C DDC port (CFG-I2COUTPROT-NMEA = 0).
 *         Use with polled UBX (neo_m9n_poll_pvt) to stop the continuous FIFO
 *         stream that maximizes I2C clock-stretch / bus-wedge exposure.
 * @return 0 on success, -1 on I2C write error.
 */
int8_t neo_m9n_disable_nmea_i2c(void);

/**
 * @brief  Poll one UBX-NAV-PVT solution (position/velocity/time) — a single
 *         bounded transaction, no continuous streaming.  Preferred over NMEA
 *         streaming on the shared bus.  Call neo_m9n_disable_nmea_i2c() once
 *         first so only UBX responses are in the FIFO.
 * @param  out  Parsed fix written here on success.
 * @return 0 on success, -1 if no valid PVT frame was read.
 */
int8_t neo_m9n_poll_pvt(neo_m9n_pvt_t *out);

/**
 * @brief  Enable periodic UBX-NAV-PVT output on the I2C port at `rate` (output
 *         every `rate` nav epochs; 1 = every epoch ≈ 1 Hz).  Use with
 *         neo_m9n_disable_nmea_i2c() + neo_m9n_read_pvt() — the robust pattern
 *         for the shared bus (steady low-volume traffic, no poll-response
 *         timing).  Call once after neo_m9n_init().
 * @return 0 on success, -1 on I2C write error.
 */
int8_t neo_m9n_enable_ubx_pvt(uint8_t rate);

/**
 * @brief  Non-blocking: drain the FIFO into an internal accumulator and, when
 *         a complete UBX-NAV-PVT frame has arrived, parse it into @p out.
 *         Call every GPS slot.  Frames span multiple calls (FIFO fills at the
 *         nav rate), so most calls return -1; a call returns 0 on the tick a
 *         fresh frame completes.  Requires neo_m9n_enable_ubx_pvt() first.
 * @return 0 if a fresh PVT was parsed this call, -1 otherwise.
 */
int8_t neo_m9n_read_pvt(neo_m9n_pvt_t *out);

/**
 * @brief  Read n bytes from the NEO-M9N stream register (0xFF).
 *         Call neo_m9n_bytes_available() first and only request as many
 *         bytes as are available; overreading returns 0xFF padding bytes.
 * @param  buf  Output buffer of at least n bytes.
 * @param  n    Number of bytes to read.
 * @return 0 on success, -1 on I2C error.
 */
int8_t  neo_m9n_read(uint8_t *buf, uint16_t n);

/**
 * @brief  Read up to max_len bytes from the stream, stopping at the first
 *         0xFF sentinel (u-blox DDC returns 0xFF when the FIFO is empty).
 *         Prefer this over neo_m9n_bytes_available() — the 0xFD counter
 *         does not update reliably on the NEO-M9N.
 * @return Number of valid bytes written to buf, or -1 on I2C error.
 */
int16_t neo_m9n_read_stream(uint8_t *buf, uint16_t max_len);

#endif /* NEO_M9N_H */
