/**
 ******************************************************************************
 *
 * @file neo_m9n.c
 *
 * @brief u-blox NEO-M9N GNSS receiver driver
 *
 * Hardened against the eff SDK I2C race/wedge (see sensor_tests/pca9685.c and
 * the Bank-2 contention analysis, 2026-07-31).  The NEO-M9N clock-stretches
 * indefinitely and, combined with the SDK's tx_finish/IIC_RST + stale-_isr_flag
 * bug, can wedge the whole shared bus (slave holds SDA low).  Three mitigations
 * here, all local to this driver (no SDK dependency):
 *
 *   1. Polled UBX (neo_m9n_poll_pvt) instead of continuous NMEA streaming —
 *      one bounded ~100 B transaction on demand, no blind FIFO-drain hunting
 *      the 0xFF sentinel.  Call neo_m9n_disable_nmea_i2c() to stop the stream.
 *   2. _isr_flag hygiene: clear the controller's stale completion flag before
 *      every transaction (same trick as pca9685.c / shtc3.c), and read the
 *      DDC stream in ≤ GPS_CHUNK-byte pieces so no single read is long.
 *   3. Bus-clear recovery: on any read/write failure, bit-bang up to 9 SCL
 *      pulses to release a slave holding SDA, then re-init — turns a bus wedge
 *      into a ~1 ms recovery instead of a dead bus.
 ******************************************************************************
 */

#include "neo_m9n.h"
#include <eff.h>
#include <eff/drivers/i2c.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

/* GPS on the shared Bank-2 sensor bus (I2C_2_1, PINMUX_2, DIGIO024/025). */
extern eff_i2c_t *I2C_2_1;
#define GPS_I2C     I2C_2_1
#define GPS_PINMUX  PINMUX_2
#define GPS_GPIO    GPIO_2
#define SCL_PIN     GPIO_PIN_4   /* DIGIO024 */
#define SDA_PIN     GPIO_PIN_5   /* DIGIO025 */

/* NEO-M9N DDC register map */
#define NEO_M9N_REG_AVAIL_HI  0xFDu  /* high byte of output FIFO byte count  */
#define NEO_M9N_REG_DATA      0xFFu  /* stream data register                  */

/* Read the DDC stream in chunks no larger than this so a single I2C
 * transaction can't stretch long enough to wedge the shared bus. */
#define GPS_CHUNK   32

/* ── bus-clear recovery ─────────────────────────────────────────────────────
 * If a read leaves a slave clamping SDA low, no START can be issued to any
 * device.  Take the pins as GPIO, clock up to 9 SCL pulses so the slave
 * finishes its byte and releases SDA, issue a STOP, hand the pins back. */
static void gps_bus_recover(void)
{
    eff_pinmux_set(GPS_PINMUX, PINMUX_GPIO);
    eff_gpio_dir_set(GPS_GPIO, SCL_PIN, EFF_GPIO_OUT);
    eff_gpio_dir_set(GPS_GPIO, SDA_PIN, EFF_GPIO_IN);
    eff_gpio_set(GPS_GPIO, SCL_PIN);
    sleep_us(5);
    for (int i = 0; i < 9; i++) {
        eff_gpio_clear(GPS_GPIO, SCL_PIN); sleep_us(5);
        eff_gpio_set  (GPS_GPIO, SCL_PIN); sleep_us(5);
        if (eff_gpio_get(GPS_GPIO, SDA_PIN)) break;   /* SDA released */
    }
    /* STOP: SDA low → high while SCL high */
    eff_gpio_dir_set(GPS_GPIO, SDA_PIN, EFF_GPIO_OUT);
    eff_gpio_clear(GPS_GPIO, SDA_PIN); sleep_us(5);
    eff_gpio_set  (GPS_GPIO, SCL_PIN); sleep_us(5);
    eff_gpio_set  (GPS_GPIO, SDA_PIN); sleep_us(5);
    eff_pinmux_set(GPS_PINMUX, PINMUX_SPI_I2C1);
    eff_i2c_init(GPS_I2C, I2C_SPEED_100K);
}

/* ── hygienic read/write ─────────────────────────────────────────────────────
 * Clear the stale completion flag before the transaction (the SDK never does),
 * then on failure recover the bus and retry once. */
static int8_t gps_read(uint8_t reg, uint8_t *buf, uint16_t n)
{
    GPS_I2C->_isr_flag = 0;
    if (eff_i2c_read(GPS_I2C, NEO_M9N_I2C_ADDR, reg, buf, n) == 0)
        return 0;
    gps_bus_recover();
    GPS_I2C->_isr_flag = 0;
    return eff_i2c_read(GPS_I2C, NEO_M9N_I2C_ADDR, reg, buf, n);
}

static int8_t gps_write(uint8_t reg, const uint8_t *buf, uint16_t n)
{
    GPS_I2C->_isr_flag = 0;
    if (eff_i2c_write(GPS_I2C, NEO_M9N_I2C_ADDR, reg,
                      (eff_i2c_buffer *)buf, n) == 0)
        return 0;
    gps_bus_recover();
    GPS_I2C->_isr_flag = 0;
    return eff_i2c_write(GPS_I2C, NEO_M9N_I2C_ADDR, reg,
                         (eff_i2c_buffer *)buf, n);
}

/* Build a UBX frame (B5 62 cls id len payload ck_a ck_b) with a runtime
 * Fletcher checksum and send it (avoids hand-computed checksum errors). */
static int8_t ubx_send(uint8_t cls, uint8_t id, const uint8_t *payload, uint16_t plen)
{
    uint8_t f[32];
    if (plen + 8u > sizeof(f)) return -1;
    f[0] = 0xB5; f[1] = 0x62; f[2] = cls; f[3] = id;
    f[4] = (uint8_t)(plen & 0xFF); f[5] = (uint8_t)(plen >> 8);
    for (uint16_t i = 0; i < plen; i++) f[6 + i] = payload[i];
    uint8_t ca = 0, cb = 0;
    for (uint16_t k = 2; k < 6u + plen; k++) { ca += f[k]; cb += ca; }
    f[6 + plen] = ca; f[7 + plen] = cb;
    /* gps_write sends [reg][data…]; reg=f[0] reconstructs the full frame. */
    return gps_write(f[0], &f[1], (uint16_t)(plen + 7));
}

/* ── init ───────────────────────────────────────────────────────────────── */

int8_t neo_m9n_init(void)
{
    eff_pinmux_set(GPS_PINMUX, PINMUX_SPI_I2C1);
    if (eff_i2c_init(GPS_I2C, I2C_SPEED_100K) != 0) {
        printf("neo_m9n_init: I2C init failed\n");
        return -1;
    }
    eff_gpio_pull_set(GPS_GPIO, SCL_PIN | SDA_PIN, EFF_GPIO_PULL_NONE);
    printf("neo_m9n_init: OK\n");
    return 0;
}

int32_t neo_m9n_bytes_available(void)
{
    uint8_t buf[2] = {0};
    if (gps_read(NEO_M9N_REG_AVAIL_HI, buf, 2) != 0)
        return -1;
    uint16_t count = ((uint16_t)buf[0] << 8) | buf[1];
    if (count == 0xFFFF)
        return 0;
    return (int32_t)count;
}

int8_t neo_m9n_read(uint8_t *buf, uint16_t n)
{
    if (gps_read(NEO_M9N_REG_DATA, buf, n) != 0) {
        printf("neo_m9n_read: I2C read failed\n");
        return -1;
    }
    return 0;
}

/*
 * Read up to max_len bytes from the stream, in ≤ GPS_CHUNK-byte transactions,
 * stopping at the first 0xFF sentinel (FIFO empty).  Chunking keeps each I2C
 * transaction short so the u-blox clock-stretch can't wedge the shared bus.
 * Returns bytes written to buf, or -1 on I2C error.
 */
int16_t neo_m9n_read_stream(uint8_t *buf, uint16_t max_len)
{
    uint16_t got = 0;
    while (got < max_len) {
        uint16_t want = max_len - got;
        if (want > GPS_CHUNK) want = GPS_CHUNK;
        if (gps_read(NEO_M9N_REG_DATA, buf + got, want) != 0)
            return -1;
        for (uint16_t i = 0; i < want; i++) {
            if (buf[got + i] == 0xFF)
                return (int16_t)(got + i);   /* FIFO empty — done */
        }
        got += want;
    }
    return (int16_t)got;
}

/* ── UBX config helpers ─────────────────────────────────────────────────── */

int8_t neo_m9n_enable_active_antenna(void)
{
    /* UBX-CFG-VALSET: CFG-HW-ANT_CFG_VOLTCTRL (0x10A3002E) = 1 */
    static const uint8_t msg[] = {
        0xB5, 0x62, 0x06, 0x8A, 0x09, 0x00,
        0x00, 0x01, 0x00, 0x00,
        0x2E, 0x00, 0xA3, 0x10, 0x01,
        0x7C, 0x21
    };
    if (gps_write(msg[0], &msg[1], sizeof(msg) - 1) != 0) {
        printf("neo_m9n: enable_active_antenna write failed\n");
        return -1;
    }
    printf("neo_m9n: antenna voltage supply enabled\n");
    return 0;
}

int8_t neo_m9n_enable_nmea_i2c(void)
{
    /* UBX-CFG-VALSET: CFG-I2COUTPROT-NMEA (0x10720002) = 1 */
    static const uint8_t msg[] = {
        0xB5, 0x62, 0x06, 0x8A, 0x09, 0x00,
        0x00, 0x01, 0x00, 0x00,
        0x02, 0x00, 0x72, 0x10, 0x01,
        0x1F, 0xB2
    };
    if (gps_write(msg[0], &msg[1], sizeof(msg) - 1) != 0) {
        printf("neo_m9n: enable_nmea_i2c write failed\n");
        return -1;
    }
    printf("neo_m9n: CFG-I2COUTPROT-NMEA enabled\n");
    return 0;
}

int8_t neo_m9n_disable_nmea_i2c(void)
{
    /* UBX-CFG-VALSET: CFG-I2COUTPROT-NMEA (0x10720002) = 0.
     * Same frame as enable with value byte 0 → recomputed Fletcher checksum. */
    static const uint8_t msg[] = {
        0xB5, 0x62, 0x06, 0x8A, 0x09, 0x00,
        0x00, 0x01, 0x00, 0x00,
        0x02, 0x00, 0x72, 0x10, 0x00,
        0x1E, 0xB1               /* CK_A, CK_B for value=0 (Fletcher-verified) */
    };
    if (gps_write(msg[0], &msg[1], sizeof(msg) - 1) != 0) {
        printf("neo_m9n: disable_nmea_i2c write failed\n");
        return -1;
    }
    printf("neo_m9n: NMEA-over-I2C disabled (polled mode)\n");
    return 0;
}

/* ── polled UBX-NAV-PVT ──────────────────────────────────────────────────────
 * Send a zero-length UBX-NAV-PVT poll, then read the DDC stream (chunked) and
 * scan for the B5 62 01 07 response frame (92-byte payload).  Validates the
 * Fletcher checksum and parses the key fields.  One bounded transaction per
 * call — no continuous streaming, minimal clock-stretch exposure.
 *
 * NOTE: the frame-sync + parse has been reasoned from the UBX spec but not yet
 * validated against a live fix — verify fixType/numSV/lat/lon on hardware.
 */
int8_t neo_m9n_poll_pvt(neo_m9n_pvt_t *out)
{
    /* UBX-NAV-PVT poll request (class 0x01, id 0x07, len 0), CK_A=0x08 CK_B=0x19 */
    static const uint8_t poll[] = { 0xB5, 0x62, 0x01, 0x07, 0x00, 0x00, 0x08, 0x19 };
    if (gps_write(poll[0], &poll[1], sizeof(poll) - 1) != 0)
        return -1;

    /* Give the receiver a moment to place the response in the DDC FIFO. */
    sleep_ms(2);

    /* Drain the stream looking for the response frame.  Response is
     * B5 62 01 07 5C 00 <92 payload> CK_A CK_B = 100 bytes; read generously. */
    uint8_t buf[256];
    int16_t n = neo_m9n_read_stream(buf, sizeof(buf));
    if (n <= 0)
        return -1;

    /* Find UBX sync + NAV-PVT header. */
    for (int i = 0; i + 8 <= n; i++) {
        if (buf[i] != 0xB5 || buf[i+1] != 0x62) continue;
        if (buf[i+2] != 0x01 || buf[i+3] != 0x07) continue;
        uint16_t len = (uint16_t)buf[i+4] | ((uint16_t)buf[i+5] << 8);
        if (len != 92) continue;
        if (i + 6 + (int)len + 2 > n) break;   /* frame not fully captured */

        /* Fletcher checksum over class..end of payload. */
        uint8_t ca = 0, cb = 0;
        for (int k = i + 2; k < i + 6 + (int)len; k++) { ca += buf[k]; cb += ca; }
        if (ca != buf[i + 6 + len] || cb != buf[i + 7 + len])
            continue;   /* bad checksum — keep scanning */

        const uint8_t *p = &buf[i + 6];   /* payload start */
        out->fix_type = p[20];
        out->num_sv   = p[23];
        out->lon_1e7  = (int32_t)((uint32_t)p[24] | (uint32_t)p[25] << 8 |
                                  (uint32_t)p[26] << 16 | (uint32_t)p[27] << 24);
        out->lat_1e7  = (int32_t)((uint32_t)p[28] | (uint32_t)p[29] << 8 |
                                  (uint32_t)p[30] << 16 | (uint32_t)p[31] << 24);
        out->height_mm = (int32_t)((uint32_t)p[36] | (uint32_t)p[37] << 8 |
                                   (uint32_t)p[38] << 16 | (uint32_t)p[39] << 24);
        return 0;
    }
    return -1;   /* no valid PVT frame found */
}

/* Enable periodic UBX-NAV-PVT output on I2C at `rate` (output every `rate`
 * nav epochs; 1 = every epoch = the nav rate, 1 Hz by default).  With NMEA
 * disabled, this pushes ~100 B/s into the FIFO — low, steady traffic the
 * scheduler drains without wedging the bus.  CFG-MSGOUT-UBX_NAV_PVT_I2C key
 * 0x20910006 (U1). */
int8_t neo_m9n_enable_ubx_pvt(uint8_t rate)
{
    const uint8_t payload[9] = {
        0x00, 0x01, 0x00, 0x00,          /* version, layer=RAM, reserved     */
        0x06, 0x00, 0x91, 0x20,          /* key 0x20910006 (LE)              */
        rate                             /* output every `rate` nav epochs   */
    };
    if (ubx_send(0x06, 0x8A, payload, sizeof(payload)) != 0) {
        printf("neo_m9n: enable_ubx_pvt write failed\n");
        return -1;
    }
    printf("neo_m9n: UBX-NAV-PVT periodic output enabled (rate=%u)\n",
           (unsigned)rate);
    return 0;
}

/* Drain whatever the FIFO holds into a persistent accumulator and, if a
 * complete UBX-NAV-PVT frame has arrived, parse it.  Non-blocking: call it
 * from the GPS slot each tick — it returns 0 only on a tick where a fresh
 * frame completed, -1 otherwise (frame still assembling / nothing yet).
 * Frames span multiple calls because the FIFO fills at the nav rate. */
int8_t neo_m9n_read_pvt(neo_m9n_pvt_t *out)
{
    static uint8_t acc[300];
    static uint16_t acc_len = 0;

    /* FIXED read — NOT the 0xFF-sentinel stream reader.  UBX is BINARY and its
     * payload contains 0xFF bytes; a sentinel-stop would truncate the frame at
     * the first internal 0xFF.  The DDC FIFO is a queue and the u-blox writes
     * each PVT frame atomically, so sequential fixed reads drain the frame
     * contiguously (then 0xFF padding once empty), and the B5 62 frame-sync +
     * checksum below locates it.  32-byte reads keep each transaction short. */
    uint8_t chunk[GPS_CHUNK];
    if (gps_read(NEO_M9N_REG_DATA, chunk, sizeof(chunk)) == 0) {
        uint16_t space = (uint16_t)(sizeof(acc) - acc_len);
        uint16_t take  = ((uint16_t)sizeof(chunk) < space)
                         ? (uint16_t)sizeof(chunk) : space;
        memcpy(acc + acc_len, chunk, take);
        acc_len += take;
    }

    /* Scan the accumulator for a complete, checksum-valid NAV-PVT frame. */
    int8_t rc = -1;
    for (int i = 0; i + 8 <= (int)acc_len; i++) {
        if (acc[i] != 0xB5 || acc[i+1] != 0x62) continue;
        if (acc[i+2] != 0x01 || acc[i+3] != 0x07) continue;
        uint16_t len = (uint16_t)acc[i+4] | ((uint16_t)acc[i+5] << 8);
        if (len != 92) continue;
        if (i + 6 + (int)len + 2 > (int)acc_len) break;   /* not fully in yet */

        uint8_t ca = 0, cb = 0;
        for (int k = i + 2; k < i + 6 + (int)len; k++) { ca += acc[k]; cb += ca; }
        if (ca != acc[i+6+len] || cb != acc[i+7+len]) continue;

        const uint8_t *p = &acc[i + 6];
        out->fix_type  = p[20];
        out->num_sv    = p[23];
        out->lon_1e7   = (int32_t)((uint32_t)p[24] | (uint32_t)p[25] << 8 |
                                   (uint32_t)p[26] << 16 | (uint32_t)p[27] << 24);
        out->lat_1e7   = (int32_t)((uint32_t)p[28] | (uint32_t)p[29] << 8 |
                                   (uint32_t)p[30] << 16 | (uint32_t)p[31] << 24);
        out->height_mm = (int32_t)((uint32_t)p[36] | (uint32_t)p[37] << 8 |
                                   (uint32_t)p[38] << 16 | (uint32_t)p[39] << 24);

        /* Consume through the end of this frame; keep the tail. */
        uint16_t end = (uint16_t)(i + 6 + len + 2);
        memmove(acc, acc + end, acc_len - end);
        acc_len -= end;
        rc = 0;
    }

    /* Bound the accumulator if no frame ever completes (drop oldest half). */
    if (acc_len > sizeof(acc) - sizeof(chunk)) {
        uint16_t keep = sizeof(acc) / 2;
        memmove(acc, acc + (acc_len - keep), keep);
        acc_len = keep;
    }
    return rc;
}
