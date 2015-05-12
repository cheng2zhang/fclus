#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 3.5 font "Times, 24"
set output "ljmix1.eps"



set xtics 5.0 font "Times,20" offset 0, 0.5
set mxtics 5
set ytics 1.0 font "Times,20" offset 0.5, 0
set mytics 2

set lmargin 5.5
set rmargin 1
set tmargin 0.7
set bmargin 2.3

color1  = "#224488"
color2  = "#aa0000"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 1.8
set style line 2  lw 2.0 lt 2 lc rgb color2  pt 8   ps 2.0


set xlabel "Cluster size" offset 0, 0.9
set ylabel "Potential of mean force" offset 1.0, 0
set key left bottom Left reverse spacing 2.0 font "Times, 24"

plot [1:27][-4.2:] \
  "../../data/ljmix/sample1/vcls.dat" u 1:(-$2) w l ls 1 t "{/Symbol-Oblique r} = 0.5", \
  "../../data/ljmix/sample2/vcls.dat" u 1:(-$2) w l ls 2 t "{/Symbol-Oblique r} = 0.4"



unset output
set terminal pop
reset
