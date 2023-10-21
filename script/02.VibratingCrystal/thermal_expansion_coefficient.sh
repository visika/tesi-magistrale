#!/usr/bin/env sh
python 3-point_derivative.py
FILE=derivative.s4.k2
paste temp+equilibrium_volume temp+equilibrium_volume.derivative | awk '{print $1,$2,$4}' >$FILE
gnuplot alpha.g
