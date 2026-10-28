load "plots/common.gp"
load "data/a1_tomo_params.gp"
dt_tomo = delta_star; lt_tomo = lambda_star
load "data/a1_regr_params.gp"
dt_regr = delta_star; lt_regr = lambda_star
set output "figures/a1_delta_lambda.svg"
set xlabel "λ"
set ylabel "δ_{min}(λ)"
set xrange [0.85:1.2]
set yrange [0:0.25]
set key top center
set label 1 sprintf("δ* = %.4f", dt_tomo) at lt_tomo, dt_tomo offset 1,-1.2
set label 2 sprintf("δ* = %.4f", dt_regr) at lt_regr, dt_regr offset 1,0.8
plot "data/a1_tomo_delta_lambda.dat" u 1:2 w l lw 2 lc rgb "#1f77b4" t "томография: rad A_1 = δ·(1 1)", \
     "data/a1_regr_delta_lambda.dat" u 1:2 w l lw 2 lc rgb "#d62728" t "регрессия: rad A_1 = δ·(1 0)", \
     "+" u (lt_tomo):(dt_tomo) every ::0::0 w p pt 7 ps 1.2 lc rgb "#1f77b4" notitle, \
     "+" u (lt_regr):(dt_regr) every ::0::0 w p pt 7 ps 1.2 lc rgb "#d62728" notitle
