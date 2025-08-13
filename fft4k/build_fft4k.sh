#!bin/bash

BLD_DIR=bld

effcc -O3 -flto expected.c -o $BLD_DIR/expected.o -c -sim
effcc -O3 -flto fft4k.c -o $BLD_DIR/fft4k.o -c -sim
effcc -O3 -flto main.c -o $BLD_DIR/main.o -c -sim
effcc -O3 -flto $BLD_DIR/main.o $BLD_DIR/fft4k.o $BLD_DIR/expected.o -o $BLD_DIR/fft4k -sim
