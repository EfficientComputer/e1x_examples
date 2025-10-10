#!/bin/bash

mkdir bld
BLD_DIR=bld

effcc -O3 -flto generator_80211.c -o $BLD_DIR/generator_80211.o -c -sim
effcc -O3 -flto parity_80211.c -o $BLD_DIR/parity_80211.o -c -sim
effcc -O3 -flto ldpc.c -o $BLD_DIR/ldpc.o -c -sim
effcc -O3 -flto main.c -o $BLD_DIR/main.o -c -sim
effcc -O3 -flto $BLD_DIR/ldpc.o $BLD_DIR/parity_80211.o $BLD_DIR/generator_80211.o $BLD_DIR/main.o -o $BLD_DIR/ldpc -sim