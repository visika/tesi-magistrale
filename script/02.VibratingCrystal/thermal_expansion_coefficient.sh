#!/usr/bin/env sh
python 3-point_derivative.py
paste temp+equilibrium_volume temp+equilibrium_volume.derivative | awk '{print $1,$2,$4}' > derivative.s4.k2
gnuplot alpha.g
