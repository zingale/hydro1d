# plot the Sod results.  Here we assume that we have a file
# called sodx.out

set term pngcairo size 600,800
set output 'sod.png'

set multiplot;

set size 1, 0.32;

set xlabel "x";

set style line 1  lw 1 lc 3

set origin 0.0, 0.666;
set ylabel "density";
plot 'sod_0070' using 1:2 notitle with points lc 1,\
     'sod-exact.out' using 1:2 notitle with lines ls 1;

set origin 0.0, 0.333;
set ylabel "velocity";
plot 'sod_0070' using 1:5 notitle with points lc 1,\
     'sod-exact.out' using 1:3 notitle with lines ls 1;

set origin 0.0, 0.0;
set ylabel "pressure";
plot 'sod_0070' using 1:6 notitle with points lc 1,\
     'sod-exact.out' using 1:4 notitle with lines ls 1;

unset multiplot;
set term x11;