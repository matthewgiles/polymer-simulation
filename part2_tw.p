set key font ",15"
set xlabel font ",15"
set ylabel font ",15"
set tics font ", 12"
set xlabel "Timesteps"
set ylabel "Tw"
set xrange[0:150000]
plot "average_dynamics_8_30" u ($1-7000000):3 w lines title "E_{Tw} = 30" lc 1, \
     "average_dynamics_8_50" u ($1-7000000):3 w lines title "E_{Tw} = 50" lc 0