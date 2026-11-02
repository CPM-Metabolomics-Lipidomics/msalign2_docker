#!/bin/bash
#${1} is the .alignment file
#${2} is the .plf file

#--- GNU Plot script starts here -------------
gnuplot << EOF
set terminal png size 1280,960
set output "${1}.png"
set title "${1}"
set xrange [0: ]
set yrange [0: ]
set xlabel "Retention time datafile 1 [s]"
set ylabel "Retention time datafile 2 [s]"
set pointsize 2
plot "${1}" using 3:6 notitle with points,\
     "${2}" using 1:2 notitle with line,\
     "${2}" using 1:2 notitle with points

EOF
