__author__ = 'troy'

from scipy import io
import scipy.sparse as sps
import numpy as np
from math import exp, log
from gurobipy import *

from generalClasses import *


# This is the class for each scenario. Essentially it is the scenario dose and overall dose for that scenario
# along with the beamlet intensities for the second stage for that scenario
class imrt_scenario(object):
    def __init__(self, data, num, m, z1dose):
        assert (isinstance(data, imrt_data))

        #Scenario index
        self.num = num
        print 'building scenario', self.num
        # build dose variables
        self.z2 = [m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z2_' + str(self.num) + '_' + str(j)) for j in
                   xrange(data.nVox)]
        m.update()
        #initialize dose constraint
        self.doseConstr2 = [m.addConstr(-self.z2[j], GRB.EQUAL, 0) for j in xrange(data.nVox)]
        m.update()
        #add in beamlet intensities
        self.x2 = [m.addVar(lb=0, vtype=GRB.CONTINUOUS,
                            column=Column(np.array(data.Dmat.getrow(i).todense()).flatten().tolist(), self.doseConstr2),
                            name='x2_' + str(self.num) + '_' + str(i))
                   for i in xrange(data.nBix)]
        m.update()
        # Add overall dose variable
        self.zS = [m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in xrange(data.nVox)]
        m.update()
        # add linking player
        self.zLinkingConstraint = [m.addConstr(self.zS[j], GRB.EQUAL,
                                               data.stageOneFraction * z1dose[j] + (1 - data.stageOneFraction) *
                                               self.z2[j], name="zlinking_" + str(self.num) + "_" + str(j)) for j in
                                   range(data.nVox)]
        m.update()


# THIS IS TROY'S ADAPTIVE LUNG CLASS. VICTOR WILL NEED TO WRITE HIS OWN OF THESE FOR LIVER
# However, this is probably a decent template, so pay attention
# The general idea is to have all of your method-specific data and variables in here. For me,
# that is the PTV EUDs and mean lung doses
class imrt_adaptiveLung(object):
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
            #already taken care of in the hard scenario lung bounds
            pass
        else:
            print 'Invalid option parameter'
            exit()

    # This method saves the objective function, bounds those values, then minimize dose to all voxels
    def cleanup(self, data, m, scenarios, structures):
        # get optimal PTV bounds
        bestPTVBounds = [self.ptvEUD[s].X for s in range(self.nscen)]

        # set optimal PTV bounds
        for s in range(self.nscen):
            self.ptvEUD[s].setAttr('LB', bestPTVBounds[s])
        m.update()

        # reset objective
        nullObj = LinExpr()
        nullObj += 0.0
        m.setObjective(nullObj)
        m.update()

        # build cleanup objective
        for s in range(self.nscen):
            for j in range(data.nVox):
                scenarios[s].zS[j].setAttr('Obj', 1)
        m.update()


        # set solver to minimize
        m.setAttr("ModelSense", -1)
        m.update()

        # uncomment these lines to output the .lp file for sanity checking
        # print 'Writing out model'
        # m.write('outCleanup.lp')
        # print 'Model writing done'

        pass

    # This builds an EUD variable for the stage one dose (see passing stage-one dose vector)
    def buildPTVEUDStageOneLB(self, data, m, s1doseVec, struct, bound):
        # build EUD variable. The struct function buildEUDBound returns several variables, so [0] just takes the first one
        self.ptvs1EUDvar = struct.buildEUDBound(s1doseVec, m, 0, data, makeIfZero=True)[0]
        # set lower bound
        self.ptvs1EUDvar.setAttr('LB', bound)
        m.update()


    # This builds the linear constraints along the points that pass through the curve we're approximating
    # via PWL constraints
    # Troy can provide more information on this, but it is better to just figure out the algebra on your own instance
    def buildPWLoption1Constraints(self, data, m, scenarios):
        self.lungMeanPWLvars = [m.addVar(lb=0, vtype=GRB.CONTINUOUS) for s in range(self.nscen)]
        m.update()
        for s in range(self.nscen):
            for x in range(len(self.option1X) - 1):
                lhs = LinExpr()
                if (self.option1X[s][x] - self.option1X[s][x + 1] == 0):
                    coef = -(self.option1X[s][x + 1] - self.option1X[s][x]) / (
                        self.option1Y[s][x + 1] - self.option1Y[s][x])
                    ub = self.option1X[s][x] - ((self.option1X[s][x + 1] - self.option1X[s][x]) / (
                        self.option1Y[s][x + 1] - self.option1Y[s][x])) * (self.option1Y[s][x])
                    lhs += self.lungMeanVars[s]
                    lhs += coef * self.lungMeanPWLvars[s]
                    m.addConstr(lhs, GRB.LESS_EQUAL, ub)
                else:
                    coef = -(self.option1Y[s][x + 1] - self.option1Y[s][x]) / (
                        self.option1X[s][x + 1] - self.option1X[s][x])
                    ub = self.option1Y[s][x] - ((self.option1Y[s][x + 1] - self.option1Y[s][x]) / (
                        self.option1X[s][x + 1] - self.option1X[s][x])) * (self.option1X[s][x])
                    lhs += coef * self.lungMeanVars[s]
                    lhs += self.lungMeanPWLvars[s]
                    m.addConstr(lhs, GRB.LESS_EQUAL, ub)
        m.update()

        # now do linking constraints
        linkingExpr = LinExpr()
        for s in range(self.nscen):
            linkingExpr += self.scenprobs[s] * self.lungMeanPWLvars[s]
        m.addConstr(linkingExpr, GRB.GREATER_EQUAL, self.alpha)
        m.update()


    # This builds the PWL objective based on EUD values
    def buildObj(self, data, m, scenarios, struct):
        assert (isinstance(struct, imrt_structure))
        # gen EUDs
        self.ptvEUD = [struct.buildEUDBound(scenarios[s].zS, m, 0, data, s, makeIfZero=True)[0] for s in
                       range(self.nscen)]
        m.update()

        # set loose lower bounds (specific to my problem)
        for s in range(self.nscen):
            self.ptvEUD[s].setAttr("LB", self.ptvLooseLower)
        m.update()
        self.ptvObjPWLvars = [m.addVar(lb=0, vtype=GRB.CONTINUOUS) for s in range(self.nscen)]
        m.update()
        # build the linear constraints that use the points along the curves we're approximating with PWL functions
        for s in range(self.nscen):
            for x in range(len(self.objX) - 1):
                lhs = LinExpr()
                if (self.objX[x] - self.objX[x + 1] == 0):
                    coef = -(self.objX[x + 1] - self.objX[x]) / (self.objY[x + 1] - self.objY[x])
                    ub = self.objX[x] - (self.objX[x + 1] - self.objX[x]) / (self.objY[x + 1] - self.objY[x]) * \
                                        self.objY[x]
                    lhs += self.ptvEUD[s]
                    lhs += coef * self.ptvObjPWLvars[s]
                    m.addConstr(lhs, GRB.LESS_EQUAL, ub)
                else:
                    coef = -(self.objY[x + 1] - self.objY[x]) / (self.objX[x + 1] - self.objX[x])
                    ub = self.objY[x] - (self.objY[x + 1] - self.objY[x]) / (self.objX[x + 1] - self.objX[x]) * \
                                        self.objX[x]
                    lhs += coef * self.ptvEUD[s]
                    lhs += self.ptvObjPWLvars[s]
                    m.addConstr(lhs, GRB.LESS_EQUAL, ub)

        self.obj = LinExpr()
        for s in range(self.nscen):
            self.obj += self.scenprobs[s] * self.ptvObjPWLvars[s]
        m.setObjective(self.obj)
        m.setAttr("ModelSense", -1)
        m.update()

    # This function gets values along the obj function curve that we'll use in our PWL constraints
    def buildPWLforObj(self, data):
        print 'building pwl objective'
        # build PWL points
        self.objX, self.objY = [], []
        currentX = self.ptvStrictLower  # todo add EUD bound to data file
        while currentX < self.ptvUpper:
            self.objX.append(currentX)
            self.objY.append(self.getObjFtn(currentX))
            currentX = currentX + self.resolution
        self.ptvLooseLower = (self.objX[1] - self.objX[0]) / (self.objY[1] - self.objY[0]) * (
            self.objX[0] * (self.objY[1] - self.objY[0]) / (self.objX[1] - self.objX[0]) - self.objY[0])
        print 'loose obj bound', self.ptvLooseLower


    # This builds points along the option 1 nonlinear constraint that we'll use in the PWL approximation
    def buildPWLforOption1Constraint(self, data):
        self.option1X, self.option1Y = [], []

        for i in range(self.nscen):
            currentX = 0
            upperLimit = -self.gamma02 / (self.gamma12 + self.gamma22 * self.biomarkers[i])
            xHolder, yHolder = [], []
            while currentX < upperLimit:
                xHolder.append(currentX)
                yHolder.append(self.getOption1Function(currentX, self.biomarkers[i]))
                currentX += self.resolution
            self.option1X.append(xHolder)
            self.option1Y.append(yHolder)


    # Based on the option, set bounds on the mean lung doses (run after making lungMeanVars).
    def buildHardScenarioLungBounds(self, data, m, scenarios, struct):
        assert (isinstance(struct, imrt_structure))
        assert (isinstance(data, imrt_data))
        self.option1MeanLungBounds = []
        for i in range(self.nscen):
            self.option1MeanLungBounds.append((-self.gamma02) / (self.gamma12 + self.gamma22 * self.biomarkers[i]))
        self.option2MeanLungBounds = []
        for i in range(self.nscen):
            self.option2MeanLungBounds.append((log((1 - self.alpha) / self.alpha) - self.gamma02) / (
            self.gamma12 + self.gamma22 * self.biomarkers[i]))

        for i in range(self.nscen):
            if self.option == 1:
                self.lungMeanVars[i].setAttr('UB', self.option1MeanLungBounds[i])
            elif self.option == 2:
                self.lungMeanVars[i].setAttr('UB', self.option2MeanLungBounds[i])
        m.update()

    #calculate option 1 function (used in getting points along constraint function curve)
    def getOption1Function(self, x, biomarker):
        return 1 / (1 + exp(self.gamma02 + (self.gamma12 + self.gamma22 * biomarker) * x))


    #calculate objective function (used in getting points along obj function curve)
    def getObjFtn(self, x):
        return pow(self.s02, exp(self.beta0 - self.beta1 *x))


# Model class - This is the main object. All other objects are contained within the stochastic model
# This class is also where the gurobi model is made and passed
# Here is where you'd generate an instance of your method-specific class

class imrt_stochastic_model(object):
    def __init__(self, inputFilename, adaptiveFilename):
        # Build data object (and read in data)
        self.data = imrt_data(inputFilename, adaptiveFilename)
        assert (isinstance(self.data, imrt_data))  # makes pycharm see the instance of data and helps with development

        # initialize gurobi model (m is what you'll add variables and constraints to)
        print 'Initializing Gurobi Model'
        self.m = Model('imrt_stoch')

        # Build stage 1 gurobi variables (x,z) and dose constraint
        print 'Building Stage 1 Gurobi Variables'
        self.z1 = [self.m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z1_' + str(j)) for j in
                   xrange(self.data.nVox)]
        self.m.update()

        # Build empty dose constraint
        print 'Building Stage-one Dose Constraint'
        self.doseConstr1 = [self.m.addConstr(-self.z1[j], GRB.EQUAL, 0) for j in xrange(self.data.nVox)]
        self.m.update()

        # Add beamlet intensities along with their coefficients in dose constraint
        print 'Populating Stage-one Dose Constraint...',
        self.x1 = [self.m.addVar(lb=0, vtype=GRB.CONTINUOUS,
                                 column=Column(np.array(self.data.Dmat.getrow(i).todense()).flatten().tolist(),
                                               self.doseConstr1), name='x1_' + str(i)) for i in xrange(self.data.nBix)]
        self.m.update()
        print 'done'

        # Build list of structure objects. Note the range assumes that STRUCTURE NUMBER STARTS AT 1!
        print 'Initializing structures'
        self.structures = [imrt_structure(self.data, s) for s in range(1, self.data.nStructs + 1)]

        # Build list of scenario structures (essentially builds x,z and z=Dx for each scenario)
        print 'Initializing scenarios'
        self.scenarios = [imrt_scenario(self.data, s, self.m, self.z1) for s in range(self.data.numscenarios)]

        # This builds the structure bounds from self.data.structBounds (see buildConstraints function in structure class)
        print 'building structure-specific constraints'
        for s in range(self.data.nStructs):
            self.structures[s].buildConstraintsAdaptive(self.data, self.m, self.z1, self.scenarios)

        # THIS IS THE IMPORTANT PART - add in your class here and remove mine. Mine, upon construction, builds the necessary variables and constraints
        # for my method. I'm going to add a few other inputs to it so I can mass-run things, but those will come later
        # Use this as a template
        print 'initializing adaptive lung class'
        self.adaLung = imrt_adaptiveLung(self.data, self.m, self.scenarios, self.structures, self.z1)

        # Uncomment to write out model
        # print 'Writing out model'
        # self.m.write('out.lp')
        # print 'Model writing done'

    # Calls solver to optimize whatever the model is currently
    def callSolver(self):
        self.m.optimize()

    # This runs my adaptive-lung-specific cleanup procedure alters the model to minimize dose to all voxels while keeping previous objectives
    def initializeCleanupADA(self):
        self.adaLung.cleanup(self.data, self.m, self.scenarios, self.structures)

    # This is the outputter that prints out the outputs to a matlabe file. NOTE: this is specific to my problem (see the ADA).
    # You will need to write your own of these .X gets the value of the gurobi variable.
    def outputVariablesADA(self):
        ptvEUDs = [self.adaLung.ptvEUD[s].X for s in range(self.data.numscenarios)]
        print 'PTV EUDs', ptvEUDs
        mld = [self.adaLung.lungMeanVars[s].X for s in range(self.data.numscenarios)]
        print 'MLD', mld
        outfilename = self.data.adaptiveFilename[:-4] + '_' + str(self.adaLung.option) + '_' + '_out.mat'

        x1 = np.array([self.x1[i].X for i in range(self.data.nBix)])
        xS = np.array(
            [np.array([self.scenarios[s].x2[i].X for i in range(self.data.nBix)]) for s in
             range(self.data.numscenarios)])
        z1 = np.array([self.z1[i].X for i in range(self.data.nVox)])
        zS = np.array(
            [np.array([self.scenarios[s].zS[i].X for i in range(self.data.nVox)]) for s in
             range(self.data.numscenarios)])
        obj = self.adaLung.obj.getValue()
        io.savemat(outfilename, {'x1': x1, 'xS': xS, 'z1': z1, 'zS': zS, 'obj': obj, 'ptvEUDs': ptvEUDs, 'mld': mld})





























