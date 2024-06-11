#!/bin/bash
#${1} is the .alignment file
#${2} is the .plf file

#--- GNU Plot script starts here -------------
gnuplot << EOF
set terminal png
set output "${1}.png"
set title "${1}"
set size 2,2
set xrange [0: ]
set yrange [0: ]
set xlabel "Scans datafile 1"
set ylabel "Scans datafile 2"
set pointsize 2
plot "${1}" using 3:6 notitle with points,\
     "${2}" using 1:2 notitle with line,\
     "${2}" using 1:2 notitle with points

EOF
