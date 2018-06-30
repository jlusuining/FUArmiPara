dataf="lumi.txt"
#set terminal jpeg large font arial size 800,600
#set output "Mass.jpg"
set terminal postscript landscape enhanced color dashed lw 0.7 "Helvetica" 10
set output "Lumi-star.eps"
#----------------------------------------------------------------
set size square 1, 1 

#set title "Fig a"
#set origin 0.0, 0.0
#set logscale x
set logscale y
set mxtics 5
#set mytics 9
#set ytics nomirror
#set format x "10^{%T}"
set format y "10^{%T}"
set xlabel "t (10^3 yr)"
set ylabel "Luminosity(Lsun)"
#set xrange [0.01:100]
#set yrange [1:1E5]
plot dataf using ($1/1000):4 title "Lumi(star)" with lines lw 1
