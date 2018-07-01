dataf="m.txt"
#set terminal jpeg large font arial size 800,600
#set output "Mass.jpg"
set terminal postscript portrait enhanced color dashed lw 0.7 "Helvetica" 10
set output "Mdot_detail[filenumber].eps"
#----------------------------------------------------------------
set size square 1, 1

#set title "Fig a"
#set origin 0.0, 0.0
#set logscale x
set logscale y
set mxtics 5
#et mytics 5
#set ytics nomirror
#set format x "10^{%T}"
set format y "10^{%T}"
set xlabel "t (yr)"
set ylabel "Mdot (Msun/yr)"
#set xrange [0.01:100]
#set yrange [1:1E5]

#set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5
#set style line 1 lt 1 lw 1 pt 7 ps 1.0 

set title "Outburst [filenumber]"

plot [[tout1]:[tout2]][] dataf using 1:2 notitle with points pt 6 ps 0.5, \
  dataf using 1:2 title "Mdot (Msun/yr)" with lines
#lt -1 pi -4 pt 6 ps 0.2 lw 1
