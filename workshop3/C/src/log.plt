set term x11 1 font "helvetica,17" linewidth 1.5 persist noraise
set title "Residual"
set xlabel "Iteration"
set logscale y

plot 'data/residual' u 1:2 notitle w l
pause 1
reread
