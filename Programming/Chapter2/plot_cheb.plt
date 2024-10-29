set terminal png size 800,600
set output 'chebyshev_plot.png'
set title "Chebyshev Interpolation of Runge Function"
set xlabel "x"
set ylabel "y"
plot 'chebyshev_exact.txt' with lines lw 2 title 'Exact', \
     'chebyshev_poly_n5.txt' with lines title 'n=5', \
     'chebyshev_poly_n10.txt' with lines title 'n=10', \
     'chebyshev_poly_n15.txt' with lines title 'n=15', \
     'chebyshev_poly_n20.txt' with lines title 'n=20'