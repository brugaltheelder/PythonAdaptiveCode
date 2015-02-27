__author__ = 'troy'

from solveFunctions import *


# params (input file)

# execute run function

datafile = 'lungmpc5ProblemData.mat'
adaptivefile = 'lungmpc5_2scen_test_2.mat'

alphas = [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925]
for a in alphas:
    print 'Running stochastic model for alpha =', a
    stochSolve(datafile, adaptivefile, alpha=a)

mldbounds = [13.1864225675219, 14.1147674461116, 14.6387435937487, 15.6642262846281, 16.7711405605526, 17.9454658494049,
             18.3766104580784, 19.6721049383816, 19.7026398309836, 20.861037748849, 21.107347590444, 21.9131633629582,
             22.370002419766, 23.5253333484721]

for m in mldbounds:
    print 'Running non-adaptive for mld bound =', m
    # nonAdaSolve(datafile,adaptivefile,m)




