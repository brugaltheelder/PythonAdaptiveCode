__author__ = 'troy'

from solveFunctions import *


# params (input file)

# execute run function

biomarkerfilename = 'fullbiomarkers.mat'

# datafile = 'lungmpc5ProblemData.mat'
# adaptivefile = 'lungmpc5_2scen_test_1.mat'
# datafile = 'lungmpc5ProblemData.mat'
# adaptivefile = 'lungmpc5_2scen_1.mat'


# Lung45 stuff
# basedir = 'lung45/'
# datafile = 'lung45ProblemData.mat'
# adaptivefile = 'lung45_4scen_equiprob_1.mat'
# adaptivefile2 = 'lung45_4scen_equiprob_2.mat'
# mldbounds = [9.18540496395532,10.1970634858365,11.6824496535941,11.8983391780155,12.8007886944296,13.1863682931844,13.7244749119663,14.5313923249729,15.142511874965,15.2642825223753,16.7069807449104,18.0908595607816,19.3866057196226,20.6400477456213]

# datafile = 'lung45ProblemData.mat'
# adaptivefile = 'lung45_4scen_uniform_1.mat'
# adaptivefile2 = 'lung45_4scen_uniform_2.mat'
# mldbounds = [9.18540496395532,10.1970634858365,11.6752436183443,11.6824496535941,12.8007886944296,12.8228752186577,13.7244749119663,14.5313923249729,14.5362200055322,15.2642825223753,15.869030247037,17.0165269358018,18.0697095244472,19.0823587078719]


#lung49 stuff
basedir = 'lung49/'
datafile = 'lung49ProblemData.mat'

# adaptivefile = 'lung49_4scen_equiprob_1.mat'
# adaptivefile2 = 'lung49_4scen_equiprob_2.mat'
# mldbounds = [9.18540496395532,10.1970634858365,11.6824496535941,11.8983391780155,12.8007886944296,13.1863682931844,13.7244749119663,14.5313923249729,15.142511874965,15.2642825223753,16.7069807449104,18.0908595607816,19.3866057196226,20.6400477456213]

adaptivefile = 'lung49_4scen_uniform_1.mat'
adaptivefile2 = 'lung49_4scen_uniform_2.mat'
mldbounds = [9.18540496395532, 10.1970634858365, 11.6752436183443, 11.6824496535941, 12.8007886944296, 12.8228752186577,
             13.7244749119663, 14.5313923249729, 14.5362200055322, 15.2642825223753, 15.869030247037, 17.0165269358018,
             18.0697095244472, 19.0823587078719]




#CHANGE SECOND ADAPTIVE FILE

alphas = [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925]
for a in alphas:
    print 'Running stochastic model option 1 for alpha =', a
    stochSolve(datafile, adaptivefile, alpha=a, caseloc=basedir)

for a in alphas:
    print 'Running stochastic model option 2 for alpha =', a
    stochSolve(datafile, adaptivefile2, alpha=a, caseloc=basedir)

# u = 0.5, 1.5
# mldbounds = [13.1864225675219, 14.1147674461116, 14.6387435937487, 15.6642262846281, 16.7711405605526, 17.9454658494049,
# 18.3766104580784, 19.6721049383816, 19.7026398309836, 20.861037748849, 21.107347590444, 21.9131633629582,
#              22.370002419766, 23.5253333484721]

# u = 0.795323131, 4.5
# mldbounds = [9.18540496395532, 10.1970634858365, 10.6752436183443, 11.6824496535941, 11.8228752186577, 12.8007886944296,
# 13.5362200055322, 13.7244749119663, 14.5313923249729, 14.869030247037, 15.2642825223753, 16.0165269358018,
#              17.0697095244472, 18.0823587078719]


for m in mldbounds:
    print 'Running non-adaptive for mld bound =', m
    nonAdaSolve(datafile, adaptivefile, m, caseloc=basedir)

exit()

import os
filenamestarts = [adaptivefile[:-4], adaptivefile2[:-4]]

for f in filenamestarts:

    dosefilenamelist = []

    for file in os.listdir("."):
        if file.endswith('out.mat') and file.startswith(f):
            dosefilenamelist.append(file)

    simulateAda(datafile, adaptivefile, dosefilenamelist, biomarkerfilename, caselocation=basedir)

filenamestarts = ['nonAda_' + f for f in filenamestarts]

for f in filenamestarts:
    dosefilenamelist = []
    for file in os.listdir("."):
        if file.endswith('out.mat') and file.startswith(f):
            dosefilenamelist.append(file)

    simulateNonAda(datafile, adaptivefile, dosefilenamelist, biomarkerfilename)
