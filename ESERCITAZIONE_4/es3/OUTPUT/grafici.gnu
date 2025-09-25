set term qt persist   # usa finestra grafica e non chiuderla subito
set multiplot layout 3,3 title "Blocchi ogni 10"

plot "pofv.dat" index 0  using 1:2 with boxes title "Blocco 0"
plot "pofv.dat" index 10 using 1:2 with boxes title "Blocco 10"
plot "pofv.dat" index 20 using 1:2 with boxes title "Blocco 20"
plot "pofv.dat" index 30 using 1:2 with boxes title "Blocco 30"
plot "pofv.dat" index 40 using 1:2 with boxes title "Blocco 40"
plot "pofv.dat" index 50 using 1:2 with boxes title "Blocco 50"
plot "pofv.dat" index 60 using 1:2 with boxes title "Blocco 60"
plot "pofv.dat" index 70 using 1:2 with boxes title "Blocco 70"
plot "pofv.dat" index 80 using 1:2 with boxes title "Blocco 80"

unset multiplot
pause -1 "Premi invio per chiudere"