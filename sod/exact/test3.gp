# plot the Sod results.  Here we assume that we have a file
# called sodx.out

#set term pngcairo size 600,800
#set output 'test3.png'

set term epscairo size 6in, 8.5in fontscale 0.75
set output 'test3.eps'

set multiplot;

set size 1, 0.32;

set xlabel "x";

set style line 1  lw 1 lc 3

set key left

set origin 0.0, 0.666;
set ylabel "density";
plot 'test3_000103' using 1:2 title 'original ppm' with points lc 1,\
     'test3_ppmT_000104' using 1:2 title 'ppm-T' with points lc 2,\
     'test3-exact.out' using 1:2 notitle with lines ls 1 lc 0;

set origin 0.0, 0.333;
set ylabel "velocity";
plot 'test3_000103' using 1:5 notitle with points lc 1,\
     'test3_ppmT_000104' using 1:5 notitle with points lc 2,\
     'test3-exact.out' using 1:3 notitle with lines ls 1 lc 0;

set origin 0.0, 0.0;
set ylabel "pressure";
plot 'test3_000103' using 1:6 notitle with points lc 1,\
     'test3_ppmT_000104' using 1:6 notitle with points lc 2,\
     'test3-exact.out' using 1:4 notitle with lines ls 1 lc 0;

unset multiplot;
set term x11;
