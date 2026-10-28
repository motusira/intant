load "plots/common.gp"
set terminal svg size 800,600 dynamic enhanced font "Liberation Serif,18" linewidth 1.5 background "white"
load sprintf("data/a1_%s_params.gp", tag)
set output sprintf("figures/a1_%s_rows.svg", tag)
set xlabel "a_{i2}"
set ylabel "a_{i1}"
set xrange [0.92:1.08]
set yrange [0.90:1.16]
set key bottom right
plot "data/a1_".tag."_rows.dat" u 1:2:($1-$3):($1+$3):($2-$4):($2+$4) w boxxyerror fs transparent solid 0.25 border lc rgb "#1f77b4" lw 2 t "строки A_1 при δ = δ*", \
     "" u 1:2 w p pt 7 ps 1.0 lc rgb "#1f77b4" t "строки mid A_1", \
     lambda_star*x w l lw 2 dt 2 lc rgb "#d62728" t sprintf("a_{i1} = %.3f·a_{i2}", lambda_star), \
     "" u 5:6 w p pt 9 ps 1.8 lc rgb "#d62728" t "строки A_1'"
