set terminal png size 800,800
set title "Heart Curve Approximated by Bézier Curves"
set xlabel "x"
set ylabel "y"
set size ratio -1
unset key

set output 'heart_plot_m10.png'
plot 'heart_m10.txt' with lines lw 2

set output 'heart_plot_m40.png'
plot 'heart_m40.txt' with lines lw 2

set output 'heart_plot_m160.png'
plot 'heart_m160.txt' with lines lw 2