set key font ",25"
set xlabel font ",25"
set ylabel font ",25"
set tics font ", 20"

set xrange [1.5:6.5]
set ylabel "<Wr> / Lk"
set xlabel "Linking Number, Lk"

plot "averages_30.dat" using 1:($6/$1):($7/$1) with errorlines title "E_{Tw} = 30" lc 1, \
     "averages_50.dat" using 1:($6/$1):($7/$1) title "E_{Tw} = 50" with errorlines lc 0