reset
set terminal epslatex color colortext standalone lw 4 size 5in, 3.5in font ',10' header \
   "\\usepackage{amsmath}\n\\usepackage{bm}\n%Set fontsize if it is defined in master document.\n\\ifx\\imglabelsize\\undefined\n\\fontsize{10}{2}\n\\else\n\\fontsize{\\imglabelsize}{2}\n\\fi"
set output 'jc_GaAs_xxxx.tex'
set xrange [0:10]
set yrange [-5:1]
set xlabel "$\\hbar \\omega \\;(\\text{eV})$"
set xtics 0.0, 2, 10 nomirror
set ytics -5.0, 2.0, 1.0 nomirror
unset key
set ylabel "$\\iota^{xxxx} (\\omega)\\; (10^{17}\\times\\text{Am/V}^3\\text{s}^2)$"
plot 'GaAs-jc_1111.dat' u ($1):(($2)/1E17) w l dt 1 lc rgb "#000000" title "$g_{\\text{jerk}}^{xxxx}$"