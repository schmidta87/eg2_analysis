set term tikz standalone size 8,6

load "palette.gnu"

set style line 1 lw 2 lc rgb colorC
set style line 2 lw 2 lc rgb colorAl
set style line 3 lw 2 lc rgb colorFe
set style line 4 lw 2 lc rgb colorPb

eventsC=9784.
eventsAl=3419.
eventsFe=11495.
eventsPb=3551.
sAl=eventsC/eventsAl
sFe=eventsC/eventsFe
sPb=eventsC/eventsPb

set yrange [0:*]
set yrange [0:800]

set key top inside right reverse Left samplen 2

set ylabel 'Counts'

sepp=0.9

set out 'leading_QSq.tex'
set xrange [1:4]
set xlabel '$Q^2$ [GeV$^2/c^2$]'
w=sepp*0.75
plot\
	"hists.txt" u 1:($2*w) index 1 w l ls 1 title 'C',\
	"hists.txt" u 1:($2*w*sAl) index 9 w l ls 2 title 'Al',\
	"hists.txt" u 1:($2*w*sFe) index 13 w l ls 3 title 'Fe$\;$',\
	"hists.txt" u 1:($2*w*sPb) index 5 w l ls 4 title 'Pb',\


set out 'leading_xB.tex'
set xrange [1.2:2]
set xlabel '$x_B$'
w=0.4*sepp
plot\
	"hists.txt" u 1:($2*w) index 3 w l ls 1 title 'C',\
	"hists.txt" u 1:($2*w*sAl) index 11 w l ls 2 title 'Al',\
	"hists.txt" u 1:($2*w*sFe) index 15 w l ls 3 title 'Fe$\;$',\
	"hists.txt" u 1:($2*w*sPb) index 7 w l ls 4 title 'Pb',\


set out 'leading_Pmiss.tex'
set xrange [0.3:1]
set xlabel '$p_{miss}$ [GeV$/c$]'
w=sepp*0.8
plot\
	"hists.txt" u 1:($2*w) index 0 w l ls 1 title 'C',\
	"hists.txt" u 1:($2*w*sAl) index 8 w l ls 2 title 'Al',\
	"hists.txt" u 1:($2*w*sFe) index 12 w l ls 3 title 'Fe$\;$',\
	"hists.txt" u 1:($2*w*sPb) index 4 w l ls 4 title 'Pb',\


set out 'leading_theta_Pmq.tex'
set xrange [100:180]
set xtics format '$%g^\circ$'
set xlabel '$\theta_{p_{miss},q}$'
w=1
plot\
	"hists.txt" u 1:($2*w) index 2 w l ls 1 title 'C',\
	"hists.txt" u 1:($2*w*sAl) index 10 w l ls 2 title 'Al',\
	"hists.txt" u 1:($2*w*sFe) index 14 w l ls 3 title 'Fe$\;$',\
	"hists.txt" u 1:($2*w*sPb) index 6 w l ls 4 title 'Pb',\


unset out
!pdflatex leading_QSq
!pdflatex leading_xB
!pdflatex leading_Pmiss
!pdflatex leading_theta_Pmq