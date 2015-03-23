__author__ = 'troy'

from scipy import io
import scipy.sparse as sps
import numpy as np
from math import exp, log
from gurobipy import *

from generalClasses import *


class imrt_model_nonada(object):
    # Initialize class
    def __init__(self, inputFilename, adaptiveFilename, mldbound=20, caselocation=''):
        self.data = imrt_data(inputFilename, adaptiveFilename, caselocation)

        assert (isinstance(self.data, imrt_data))  # makes pycharm see the instance of data and helps with development

        # initialize gurobi model (m is what you'll add variables and constraints to)
        print 'Initializing Gurobi Model'
        self.m = Model('imrt_nonAda_stoch')

        # Build z variables
        print 'Building Gurobi Variables'
        self.z = [self.m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z_' + str(j)) for j in
                  xrange(self.data.nVox)]
        self.m.update()

        # Build empty dose constraint
        print 'Building Stage-one Dose Constraint'
        self.doseConstr = [self.m.addConstr(-self.z[j], GRB.EQUAL, 0) for j in xrange(self.data.nVox)]
        self.m.update()

        # Add beamlet intensities along with their coefficients in dose constraint
        print 'Populating Stage-one Dose Constraint...',
        self.x = [self.m.addVar(lb=0, vtype=GRB.CONTINUOUS,
                                column=Column(np.array(self.data.Dmat.getrow(i).todense()).flatten().tolist(),
                                              self.doseConstr), name='x_' + str(i)) for i in xrange(self.data.nBix)]
        self.m.update()
        print 'done'

        # Build list of structure objects. Note the range assumes that STRUCTURE NUMBER STARTS AT 1!
        print 'Initializing structures'
        self.structures = [imrt_structure(self.data, s) for s in range(1, self.data.nStructs + 1)]

        # This builds the structure bounds from self.data.structBounds (see buildConstraints function in structure class)
        print 'building structure-specific constraints'
        for s in range(self.data.nStructs):
            self.structures[s].buildConstraintsNonAdaptive(self.data, self.m, self.z)

        #Build up non-adaptive class with lung bounds and objective
        print 'initializing non-adaptive lung class'
        self.nonadaLung = imrt_nonAdaptiveLung(self.data, self.m, self.structures, self.z, mldBoundValue=mldbound)

        print 'Writing initial non-adaptive out model'
        #self.m.write('outNonAda.lp')
        print 'Model writing done'


    # Calls solver to optimize whatever the model is currently
    def callSolver(self):
        self.m.optimize()

    # This runs my adaptive-lung-specific cleanup procedure alters the model to minimize dose to all voxels while keeping previous objectives
    def initializeCleanupNonADA(self):
        self.nonadaLung.cleanup(self.data, self.m, self.z, self.structures)

    # This is the outputter that prints out the outputs to a matlabe file. NOTE: this is specific to my problem (see the ADA).
    # You will need to write your own of these .X gets the value of the gurobi variable.
    def outputVariablesnonADA(self):
        ptvEUD = self.nonadaLung.ptvEUD.X
        print 'PTV EUD', ptvEUD
        mld = self.nonadaLung.lungMeanVar.X
        print 'MLD', mld
        outfilename = self.data.basedir + 'nonAda_' + self.data.adaptiveFilename[:-4] + '_' + str(
            len(self.nonadaLung.biomarkers)) + '_' + str(round(self.nonadaLung.meanLungBoundNonAda * 100)) + '_out.mat'

        x = np.array([self.x[i].X for i in range(self.data.nBix)])

        z = np.array([self.z[i].X for i in range(self.data.nVox)])

        obj = self.nonadaLung.obj.getValue()
        io.savemat(outfilename, {'x': x, 'z': z, 'obj': obj, 'ptvEUD': ptvEUD, 'mld': mld,
                                 'mldbound': self.nonadaLung.meanLungBoundNonAda})







class imrt_nonAdaptiveLung(object):
    # initizlize class
    def __init__(self, data, m, structures, zdose, mldBoundValue=20):
        # open  up data file, read into this specific structure
        # option 1 is expected bound
        # option 2 is any scenario bound
        # assert (isinstance(data, imrt_data)) # This assertion just makes pycharm work better with data, feel free to  uncomment it

        # Sets the solver to primal simplex
        m.setParam('Method', 0)

        self.meanLungBoundNonAda = mldBoundValue

        print 'Reading in adaptive lung class data'
        # Reads in problem-specific parameters
        adaMatFile = io.loadmat(data.basedir + data.adaptiveFilename)
        self.nscen = int(adaMatFile['nscen'])
        self.beta0 = float(adaMatFile['beta0'])
        self.beta1 = float(adaMatFile['beta1'])
        self.gamma02 = float(adaMatFile['gamma02'])
        self.gamma12 = float(adaMatFile['gamma12'])
        self.gamma22 = float(adaMatFile['gamma22'])
        self.s02 = float(adaMatFile['s02'])

        # update scenario's data variables
        self.biomarkers = np.array(adaMatFile['biomarkers']).flatten()
        self.scenprobs = np.array(adaMatFile['scenprob']).flatten()

        # read in which structures go to which objective
        self.ptvstruct = int(adaMatFile['ptvStruct']) - 1  # This is because of the 1 indexing in matlab
        self.lungstruct = int(adaMatFile['lungStruct']) - 1

        self.ptvStrictLower = float(adaMatFile['ptvStrictLower'])
        self.ptvLower = float(adaMatFile['ptvLower'])
        self.ptvUpper = float(adaMatFile['ptvUpper'])

        # alpha is 1- probability of RILT
        self.alpha = float(adaMatFile['alpha'])
        # option 1 is expected bound
        # option 2 is any scenario bound
        self.option = int(adaMatFile['option'])
        self.resolution = float(adaMatFile['resolution'])

        print 'Finished reading in adaptive lung data'

        ##build objective, add to model
        print 'Building Objective'
        self.buildObj(data, m, zdose, structures[self.ptvstruct])

        # build lung means
        print 'Building Lung Mean Vars'
        self.lungMeanVar = structures[self.lungstruct].buildMeanBound(zdose, m, 0)

        # Sets worst-case lung bounds
        print 'Setting hard lung bounds'
        self.buildHardLungBounds(m)


    # This method saves the objective function, bounds those values, then minimize dose to all voxels
    def cleanup(self, data, m, zDose, structures):
        # get optimal PTV bounds
        bestPTVBound = self.ptvEUD.X

        # set optimal PTV bounds

        self.ptvEUD.setAttr('LB', bestPTVBound)
        m.update()

        # reset objective
        nullObj = LinExpr()
        nullObj += 0.0
        m.setObjective(nullObj)
        m.update()

        # build cleanup objective

        for j in range(data.nVox):
            zDose[j].setAttr('Obj', 1)
        m.update()


        # set solver to minimize
        m.setAttr("ModelSense", -1)
        m.update()

        # uncomment these lines to output the .lp file for sanity checking
        print 'Writing out model'
        # m.write('outCleanupNonAda.lp')
        print 'Model writing done'

    # This builds the PWL objective based on EUD values
    def buildObj(self, data, m, zDose, struct):
        assert (isinstance(struct, imrt_structure))
        # gen EUD, set objective function to this
        self.ptvEUD = struct.buildEUDBound(zDose, m, 0, data, makeIfZero=True)[0]
        m.update()
        self.obj = LinExpr()
        self.obj += self.ptvEUD
        m.setObjective(self.obj)
        m.setAttr("ModelSense", -1)
        m.update()

    # Based on the option, set bounds on the mean lung doses (run after making lungMeanVars).
    def buildHardLungBounds(self, m):
        self.lungMeanVar.setAttr('UB', self.meanLungBoundNonAda)
        m.update()
