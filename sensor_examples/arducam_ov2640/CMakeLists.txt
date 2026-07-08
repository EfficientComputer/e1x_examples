add_eff_app(
  arducam_id
  ARCHS e1x
  TYPE
  exe
  TARGETS
  scalar
  SOURCE
  arducam_id.c)

add_eff_app(
  cam_probe
  ARCHS e1x
  TYPE
  exe
  TARGETS
  scalar
  SOURCE
  cam_probe.c)

add_eff_app(
  spi_loopback
  ARCHS e1x
  TYPE
  exe
  TARGETS
  scalar
  SOURCE
  spi_loopback.c)

add_eff_app(
  fifo_bringup
  ARCHS e1x
  TYPE
  exe
  TARGETS
  scalar
  SOURCE
  fifo_bringup.c)

add_eff_app(
  single_frame
  ARCHS e1x
  TYPE
  exe
  TARGETS
  scalar
  SOURCE
  single_frame.c)

add_eff_app(
  single_frame_jpeg_verify
  ARCHS e1x
  TYPE
  exe
  TARGETS
  scalar
  SOURCE
  single_frame_jpeg_verify.c)
