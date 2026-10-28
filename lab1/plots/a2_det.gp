load "plots/common.gp"
load "data/a2_params.gp"
set output "figures/a2_det.svg"
set xlabel "δ"
set ylabel "det"
set key top left
set xrange [0:0.25]
set arrow 1 from delta_star, graph 0 to delta_star, graph 1 nohead dt 3 lc rgb "#555555"
set label 1 sprintf("δ* = %.3f", delta_star) at delta_star, graph 0.92 offset 0.7,0
plot "data/a2_det_range.dat" u 1:2:3 w filledcurves fs transparent solid 0.2 lc rgb "#1f77b4" t "{/:Bold det} A_2(δ) = [min, max] по вершинам", \
     "" u 1:2 w l lw 2 lc rgb "#1f77b4" notitle, \
     "" u 1:3 w l lw 2 lc rgb "#1f77b4" notitle, \
     "data/a2_det_path.dat" u 1:2 w l lw 2 dt 2 lc rgb "#d62728" t "det(mid A_2 − t·y z^T)", \
     0 w l lw 1 lc rgb "black" notitle
