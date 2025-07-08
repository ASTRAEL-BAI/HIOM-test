This project is based on the original code from [hvdmaas/HIOM](https://github.com/hvdmaas/HIOM).  
The original author is hvdmaas. The code is licensed under the [GNU General Public License v3.0 (GPL-3.0)](https://www.gnu.org/licenses/gpl-3.0.html).


----------------------------------------------------------------------------------------------------------------------------------------------------------------
# HIOM
R code for "The polarization within and across individuals: the hierarchical Ising opinion model"
submitted to Journal of Complex Networks

Polarization of opinions is a societal threat. It involves psychological processes as well as group dynamics, a popular topic in statistical physics. However, the interaction between the within individual dynamics of attitude formation and across person polarization is rarely studied. By modelling individual attitudes as Ising networks of attitude elements, and approximating this behaviour by the cusp singularity, we developed a fundamentally new model of social dynamics. In this hierarchical model agents behave either discretely or continuously depending on their attention to the issue. At the individual level the model reproduces the mere thought effect and resistance to persuasion. At the social level the model implies polarization and the persuasion paradox. We propose a new intervention for escaping polarization in bounded confidence models of opinion dynamics. 

![The Black Pete simulation](figures/Anim_3.gif)

How to run
----------
Run code in Rstudio
Install packages specified in the first lines of HIOM_SW.r
add a directory figures in the directory of HIOM_SW.r

Main simulation is HIOM_SW.r (figures 4,5,6,7,8)

set pdfplot=F; PNG=F; scenario=3 to get the figure above.

set pdfplot=T; PNG=F; scenario=2 to get figure 4.

set pdfplot=T; PNG=F; scenario=3 to get figure 5.

set pdfplot=F; PNG=F; scenario=31 to get figure 6.

set pdfplot=T; PNG=F; scenario=4 to get figure 7.

set pdfplot=F; PNG=F; scenario=41 to get figure 8.



figures2&3.r produces figure 2 and 3

figureA1.r produces figure in Appendix A

----------------------------------------------------------------------------------------------------------------------------------------------------------------New:

scenario 6 : see how attention move in different dA
scenario 91 : see depolarization