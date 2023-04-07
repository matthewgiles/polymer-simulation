set key font ",15"
set xlabel font ",15"
set ylabel font ",15"
set tics font ", 10"

set xrange [25:55]
set ylabel "Variance in Tw"
set xlabel "E_{Tw}"

plot "part1_6_variance" u 1:2 title "Lk = 6, no nic" w lines
#plot "part2_8_variance" u 1:2 title "Lk = 8, with nic" w lines