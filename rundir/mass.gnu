dataf="m.txt"
#set terminal jpeg large font arial size 800,600
#set output "Mass.jpg"
set terminal postscript landscape enhanced color dashed lw 0.7 "Helvetica" 10
set output "Mass.eps"
#----------------------------------------------------------------
set size square 1, 1

#set title "Fig a"
#set origin 0.0, 0.0
#set logscale x
#set logscale y
set mxtics 5
set mytics 5
#set ytics nomirror
#set format x "10^{%T}"
#set format y "10^{%T}"
set xlabel "t (10^3 yr)"
set ylabel "Mass (Msun)"
#set xrange [0.01:100]
#set yrange [1:1E5]
plot dataf using ($1/1000):3 title "Mdisk" with lines lw 1, \
     dataf using ($1/1000):4 title "Mstar" with lines lw 1

