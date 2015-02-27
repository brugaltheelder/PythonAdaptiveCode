__author__ = 'troy'
from adaptiveClasses import *
from nonadaptiveClasses import *


def stochSolve(datafilename, adaptivefilename, alpha=-1):
    # initizlize instance of stochastic model this basically does everything other than solve the model.
    # If you want it to only add bits, then just change how the constructor works (i.e. add in
    # functions for tasks like adding the constraints, making dose, etc outside of the constructor and
    # call them from here

    stoch_mod = imrt_stochastic_model('lungmpc5ProblemData.mat', 'lungmpc5_2scen_test_1.mat', manualAlpha=alpha)

    # This is another instance that runs that is larger
    # stoch_mod = imrt_stochastic_model('lung45ProblemData.mat', 'lung45_4scen_equiprob.mat')


    # This solves whatever the model state is
    stoch_mod.callSolver()

    # This is a function that calls my cleanup subroutines for my specific class (see ADA on end?)
    stoch_mod.initializeCleanupADA()

    # Solve the new model
    stoch_mod.callSolver()

    # Output the model
    stoch_mod.outputVariablesADA()


def nonAdaSolve(datafilename, adaptivefilename, meanbound):
    nonada_mod = imrt_model_nonada(datafilename, adaptivefilename, mldbound=meanbound)
    nonada_mod.callSolver()
    nonada_mod.initializeCleanupNonADA()
    nonada_mod.callSolver()
    nonada_mod.outputVariablesnonADA()
