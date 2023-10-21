set terminal png
set output 'alpha.png'
set yrange [-10:120]
set key bottom
set title "Volumetric thermal expansion coefficient of Al"
set xlabel "Temperature (K)"
set ylabel "α · 10^6 (m^3/m^3K)"
plot 'derivative.s4.k2' u 1:($3/$2) t "VASP v.s4.k2", \
     '' u 1:($3/$2) w l ls 1 not, \
     'wilson' u ($1+273.15):($3*3) ls 2 t "Wilson's experiment", \
     '' u ($1+273.15):($3*3) w l ls 2 not
