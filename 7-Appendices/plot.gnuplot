set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output "a_n-q10.png"
set title "a_n vs Energy (E)"
set autoscale
set ylabel "a_n"
set xlabel "Energy (E)"
set key bottom left
set yrange [1.37:1.48]

plot 'anq10.dat' using 1:2:(0.5) with linespoints smooth sbezier, 'anq10.dat' using 1:2
