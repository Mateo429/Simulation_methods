set title "Evolución temporal"
set xrange[0:50]
set yrange[0:0.5]
set xlabel "tiempo"
set ylabel "Fracción poblacional"
set grid



plot "data.txt" u 1:2 w l, "data.txt" u 1:3 w l,  "data.txt" u 1:4 w l

set output "grafica.png"
replot