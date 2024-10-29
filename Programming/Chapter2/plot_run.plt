set terminal png size 800,600
set output 'runge_plot.png'
set title "Runge Function and Interpolations"
set xlabel "x"
set ylabel "y"
plot 'runge_exact.txt' with lines lw 2 title 'Exact', \
     'runge_poly_n2.txt' with lines title 'n=2', \
     'runge_poly_n4.txt' with lines title 'n=4', \
     'runge_poly_n6.txt' with lines title 'n=6', \
     'runge_poly_n8.txt' with lines title 'n=8'