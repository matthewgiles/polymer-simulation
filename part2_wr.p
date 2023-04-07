set key font ",15"
set xlabel font ",15"
set ylabel font ",15"
set tics font ", 12"
set xlabel "Timesteps"
set ylabel "Wr"
set xrange[0:1500000]
plot "average_dynamics_8_30" u ($1-7000000):4 w lines title "E_{Tw} = 30" lc 1, \
     "average_dynamics_8_50" u ($1-7000000):4 w lines title "E_{Tw} = 50" lc 0