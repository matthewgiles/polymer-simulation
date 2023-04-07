set key font ",15"
set xlabel font ",15"
set ylabel font ",15"
set tics font ", 12"
set xrange [0:35000]
set xlabel "Timesteps"
set ylabel "log(Tw)"
f(x) = a*x + b
plot "average_dynamics_11_50" u ($1-7000000):(log($3)) w lines title "E_{Tw} = 30"
fit [0:35000] f(x) "average_dynamics_11_50" u ($1-7000000):(log($3)) via a,b