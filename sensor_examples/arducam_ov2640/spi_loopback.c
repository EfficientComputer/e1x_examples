//===-- spi_loopback.c — E1x bank5 SPI loopback (jumper DIGIO052 MOSI <-> 053 MISO) ===//

#include <eff.h>
#include <eff/drivers/spi.h>

#include <stdio.h>

// SPI: bank 5 (DIGIO050 CLK, 051 CS, 052 MOSI, 053 MISO), SW10 ON.
#define LB_SPI     SPI_5
#define LB_PINMUX  PINMUX_5

int main() {
    printf("\r\nspi loopback bank5 (jumper DIGIO052 MOSI <-> DIGIO053 MISO)\r\n\r\n");

    eff_pinmux_set(LB_PINMUX, PINMUX_SPI);
    eff_spi_cfg_t c = EFF_SPI_DEFAULTS;
    c.mode = SPI_MODE_0;
    c.xfer_mode = SPI_XFER_BIDIRECTIONAL;
    c.bus_size = SPI_BUS_SINGLE;
    c.clk_div = 160;
    eff_spi_init(LB_SPI, &c);

    uint8_t pats[4] = {0x55, 0xAA, 0xA5, 0x5A};
    for (int i = 0; i < 4; i++) {
        uint8_t tx[1] = {pats[i]};
        volatile uint8_t rx[1] = {0};
        eff_spi_xfer(LB_SPI, 0, 0, tx, 1, (uint8_t *)rx, 1);
        printf("tx 0x%02x -> rx 0x%02x  %s\r\n", pats[i], rx[0],
               rx[0] == pats[i] ? "OK" : "");
    }

    return 0;
}
