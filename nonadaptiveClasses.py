__author__ = 'troy'

from scipy import io
import scipy.sparse as sps
import numpy as np
from math import exp, log
from gurobipy import *

from generalClasses import *


class imrt_model(object):
    pass


class imrt_nonAdaptiveLung(object):
    # initizlize class
    def __init__(self, data, m, scenarios, structures, z1dose):
        # open  up data file, read into this specific structure
        # option 1 is expected bound
        # option 2 is any scenario bound
        # assert (isinstance(data, imrt_data)) # This assertion just makes pycharm work better with data, feel free to  uncomment it

        # Sets the solver to primal simplex
        m.setParam('Method', 0)

        print 'Reading in adaptive lung class data'
        # Reads in problem-specific parameters
        adaMatFile = io.loadmat(data.adaptiveFilename)
        self.nscen = int(adaMatFile['nscen'])
        self.beta0 = float(adaMatFile['beta0'])
        self.beta1 = float(adaMatFile['beta1'])
        self.gamma02 = float(adaMatFile['gamma02'])
        self.gamma12 = float(adaMatFile['gamma12'])
        self.gamma22 = float(adaMatFile['gamma22'])
        self.s02 = float(adaMatFile['s02'])

        # Reads in stage-one fraction (in case you use it)
        self.stageOneFraction = float(adaMatFile['s1frac'])

        # update scenario's data variables
        self.biomarkers = np.array(adaMatFile['biomarkers']).flatten()
        self.scenprobs = np.array(adaMatFile['scenprob']).flatten()

        # read in which structures go to which objective
        self.ptvstruct = int(adaMatFile['ptvStruct']) - 1  # This is because of the 1 indexing in matlab
        self.lungstruct = int(adaMatFile['lungStruct']) - 1
        self.ptvEUDs1bound = float(adaMatFile['ptvEUDs1bound'])

        self.ptvStrictLower = float(adaMatFile['ptvStrictLower'])
        self.ptvLower = float(adaMatFile['ptvLower'])
        self.ptvUpper = float(adaMatFile['ptvUpper'])
        self.alpha = float(adaMatFile['alpha'])
        self.option = int(adaMatFile['option'])
        self.resolution = float(adaMatFile['resolution'])

        print 'Finished reading in adaptive lung data'

        ##build PWL values (points along PWL curve) for objective
        print 'Building PWL objective helper values'
        self.buildPWLforObj(data)

        # Builds points along PWL bound for option 1 (expected RILT bound)
        print 'Building PWL constraint helper values'
        self.buildPWLforOption1Constraint(data)

        ##build objective, add to model
        print 'Building Objective'
        self.buildObj(data, m, scenarios, structures[self.ptvstruct])

        # Sets a lower bound on the stage 1 ptv dose
        print 'Setting stage one bound on PTV EUD'
        self.buildPTVEUDStageOneLB(data, m, z1dose, structures[self.ptvstruct], self.ptvEUDs1bound)

        # build lung means
        print 'Building Lung Mean Vars'
        self.lungMeanVars = [structures[self.lungstruct].buildMeanBound(scenarios[s].zS, m, 0, s) for s in
                             range(self.nscen)]

        # Sets worst-case lung bounds
        print 'Setting hard lung bounds'
        self.buildHardScenarioLungBounds(data, m, scenarios, structures[self.lungstruct])

        # do the MLD bounding based on option
        print 'Building remaining lung constraints'
        if self.option == 1:
            self.buildPWLoption1Constraints(data, m, scenarios)
        elif self.option == 2:
            # already taken care of in the hard scenario lung bounds
            pass
        else:
            print 'Invalid option parameter'
            exit()