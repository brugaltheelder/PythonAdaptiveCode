__author__ = 'troy'

from scipy import io
import scipy.sparse as sps
import numpy as np
from math import exp, log
from gurobipy import *



# Data class
class imrt_data(object):
    def __init__(self, inputFilename, adaptiveFilename):
        matFile = io.loadmat(inputFilename)
        self.nVox, self.nBix, self.nDijs, self.nStructs = int(matFile['nvox']), int(matFile['nbixel']), int(
            matFile['numdijs']), int(matFile['numstructs'])
        self.numoars, self.numtargets, self.oars, self.targets = int(matFile['numoars']), int(
            matFile['numtargets']), np.array(matFile['oars']).flatten(), np.array(matFile['targets']).flatten()
        bixe = np.array(matFile['bixe2_new']).flatten() - 1
        voxe = np.array(matFile['voxe2_new_nvox']).flatten() - 1
        dijs = np.array(matFile['dijs2_new']).flatten()
        self.Dmat = sps.csr_matrix((dijs, (bixe, voxe)), shape=(self.nBix, self.nVox))
        self.maskValue = np.array(matFile['maskValue']).flatten()
        self.structPerVoxel = np.array(matFile['structs']).flatten()
        self.pickstructs = map(str, matFile['pickstructs'])

        self.stageOneFraction = 0.2  # todo read this in from data file

        # todo update these in mat file with structParams .txt file (eud weights too)
        self.structBounds = np.array(matFile['structurebounds'])
        self.structGamma = np.array(matFile['eudweights']).flatten()

        # todo read in scenario data file location and data file
        # todo add in switch for adaptive type
        self.adaptiveFilename = adaptiveFilename
        adaptiveFile = io.loadmat(self.adaptiveFilename)
        self.numscenarios = adaptiveFile['nscen']



class imrt_scenario(object):
    def __init__(self, data, num, m, z1dose):
        assert (isinstance(data, imrt_data))
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
        self.zS = [m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in xrange(data.nVox)]
        m.update()
        self.zLinkingConstraint = [m.addConstr(self.zS[j], GRB.EQUAL,
                                               data.stageOneFraction * z1dose[j] + (1 - data.stageOneFraction) *
                                               self.z2[j], name="zlinking_" + str(self.num) + "_" + str(j)) for j in
                                   range(data.nVox)]
        m.update()


class imrt_adaptiveLung(object):
    # initizlize class
    def __init__(self, data, m, scenarios, structures):
        # open  up data file, read into this specific structure
        # option 1 is expect bound
        #option 2 is any scenario bound
        assert (isinstance(data, imrt_data))

        adaMatFile = io.loadmat(data.adaptiveFilename)
        self.nscen = adaMatFile['nscen']
        self.beta0 = adaMatFile['beta0']
        self.beta1 = adaMatFile['beta1']
        self.gamma02 = adaMatFile['gamma02']
        self.gamma12 = adaMatFile['gamma12']
        self.gamma22 = adaMatFile['gamma22']
        self.s02 = adaMatFile['s02']

        #todo update matlab file with biomarker values from scenarios.txt
        # update scenario's data variables
        self.biomarkers = np.array(adaMatFile['biomarkers']).flatten()
        self.scenprobs = np.array(adaMatFile['scenprob']).flatten()

        # read in which structures go to which objective
        self.ptvstruct = adaMatFile['ptvStruct']
        self.lungstruct = adaMatFile['lungStruct']

        # todo add these to matlab file
        self.ptvStrictLower = adaMatFile['ptvStrictLower']
        self.ptvLower = adaMatFile['ptvLower']
        self.ptvUpper = adaMatFile['ptvUpper']
        self.alpha = adaMatFile['alpha']
        self.option = adaMatFile['option']

        ##build PWL values
        self.buildPWLforObj(data)
        self.buildPWLforOption1Constraint(data)

        ##build objective, add to model
        self.buildObj(data, m, scenarios, structures[self.ptvstruct])

        ##build lung means
        self.lungMeanVars = [structures[self.lungstruct].buildMeanBound(scenarios[s].zS, m, 0, s) for s in
                             range(self.nscen)]

        self.buildHardScenarioLungBounds(data, m, scenarios, structures[self.lungstruct])

        ##do the bounding based on option
        if self.option == 1:
            self.buildPWLoption1Constraints(data, m, scenarios)
        elif self.option == 2:
            pass
        else:
            print 'Invalid option parameter'
            exit()

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


    def buildObj(self, data, m, scenarios, struct):
        assert (isinstance(struct, imrt_structure))
        # gen EUDs
        self.ptvEUD = [struct.buildEUDBound(scenarios[s].zS, m, 0, data, s) for s in range(self.nscen)]
        m.update()
        # set loose lower bounds
        for s in range(self.nscen):
            self.ptvEUD[s].setAttr("LB", self.ptvLooseLower)
        m.update()
        self.ptvObjPWLvars = [m.addVar(lb=0, vtype=GRB.CONTINUOUS) for s in range(self.nscen)]
        m.update()
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
            self.obj += self.scenprobs * self.ptvObjPWLvars[s]
        m.setObjective(self.obj)
        m.setAttr("ModelSense", -1)
        m.update()

    def buildPWLforObj(self, data):

        print 'building pwl objective'
        # build PWL points
        resolution = 0.5
        self.objX, self.objY = [], []
        currentX = self.ptvStrictLower  # todo add EUD bound to data file
        while currentX < self.ptvUpper:
            self.objX.append(currentX)
            self.objY.append(self.getObjFtn(currentX))
            currentX = currentX + resolution
        self.ptvLooseLower = (self.objX[1] - self.objX[0]) / (self.objY[1] - self.objY[0]) * (
            self.objX[0] * (self.objY[1] - self.objY[0]) / (self.objX[1] - self.objX[0]) - self.objY[0])
        print 'loose obj bound', self.ptvLooseLower

    def buildPWLforOption1Constraint(self, data):
        self.option1X, self.option1Y = [], []

        for i in range(self.nscen):
            currentX = 0
            upperLimit = -self.gamma02 / (self.gamma12 + self.gamma22 * self.biomarkers[i])
            xHolder, yHolder = [], []
            while currentX < upperLimit:
                xHolder.append(currentX)
                yHolder.append(self.getOption1Function(currentX, self.biomarkers[i]))
            self.option1X.append(xHolder)
            self.option1Y.append(yHolder)


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


    def getOption1Function(self, x, biomarker):
        return 1 / (1 + exp(self.gamma02 + (self.gamma12 + self.gamma22 * biomarker) * x))


    # build objective todo write objective function generator
    def buildObjective(self, data, m, scenarios, struct):
        self.buildPWLforObj(data, struct)
        #build gurobi stuff

    def getObjFtn(self, x):
        return pow(self.s02, exp(self.beta0 - self.beta1 *x))

    # build constraints todo write constraint function generator
    def buildAdaptiveConstraint(self, data, m, scenarios, struct):
        pass


# adaLung class

# todo reads in separate data file




# Model class
class imrt_stochastic_model(object):
    def __init__(self, inputFilename):
        # Build data
        self.data = imrt_data(inputFilename)
        assert (isinstance(self.data, imrt_data))  # makes pycharm see the instance of data and helps with development
        # initialize gurobi model
        print 'Initializing Gurobi Model'
        self.m = Model('imrt_stoch')

        # Build stage 1 gurobi variables (x,z) and dose constraint
        print 'Building Stage 1 Gurobi Variables'
        self.z1 = [self.m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z1_' + str(j)) for j in
                   xrange(self.data.nVox)]
        self.m.update()
        print 'Building Stage-one Dose Constraint'
        self.doseConstr1 = [self.m.addConstr(-self.z1[j], GRB.EQUAL, 0) for j in xrange(self.data.nVox)]
        self.m.update()
        print 'Populating Stage-one Dose Constraint'
        self.x1 = [self.m.addVar(lb=0, vtype=GRB.CONTINUOUS,
                                 column=Column(np.array(self.data.Dmat.getrow(i).todense()).flatten().tolist(),
                                               self.doseConstr1), name='x1_' + str(i)) for i in xrange(self.data.nBix)]
        self.m.update()
        print 'Stage-one Dose Constraint Built'

        print 'initializing structures'
        self.structures = [imrt_structure(self.data, s) for s in range(1, self.data.nStructs + 1)]

        print 'initializing scenarios'
        self.scenarios = [imrt_scenario(self.data, s, self.m, self.z1) for s in range(self.data.numscenarios)]

        print 'building structure-specific constraints'
        for s in range(self.data.nStructs):
            self.structures[s].buildConstraints(self.data, self.m, self.z1, self.scenarios)


            #todo initizlize adaptive class

        # Uncomment to write out model
        # print 'Writing out model'
        # self.m.write('out.lp')
        # print 'Model writing done'


    def callSolver(self):
        self.m.optimize()


# Structure class
class imrt_structure(object):
    def __init__(self, data, index):
        assert (isinstance(data, imrt_data))
        self.name = data.pickstructs[
            index - 1]  # ASSUMES ANAT INDEXING STARTS AT 1 TODO FIX THIS SO IT STARTS AT 0, also below and elsewhere
        self.index = index
        self.voxels = np.where(data.structPerVoxel == index)[0]
        self.size = self.voxels.size
        self.z1bounds = data.structBounds[index - 1, 0:4]
        self.z2bounds = data.structBounds[index - 1, 4:8]
        self.zSbounds = data.structBounds[index - 1, 8:12]

    def buildConstraints(self, data, m, z1dose, scenarios):
        # for each set of bounds, for each bound value (if >0), build constraint

        print 'Generating bounds for structure number', self.index, '(', self.name, ')for z1'
        #z1bounds
        for b in range(len(self.z1bounds)):
            if self.z1bounds[b] > 0 and b == 0:
                #min constraint
                self.buildMinBound(z1dose, m, self.z1bounds[b])

            elif self.z1bounds[b] != 0 and b == 1:
                # mean constraint
                self.buildMeanBound(z1dose, m, self.z1bounds[b])

            elif self.z1bounds[b] > 0 and b == 2:
                # max constraint
                self.buildMaxBound(z1dose, m, self.z1bounds[b])

            elif self.z1bounds[b] > 0 and b == 3:
                self.z1eud = self.buildEUDBound(z1dose, m, self.z1bounds[b], data)[0]
        print 'Generating bounds for structure number', self.index, '(', self.name, ')for z2'
        # z2bounds
        for b in range(len(self.z2bounds)):
            if self.z2bounds[b] > 0 and b == 0:
                #min constraint
                for s in range(data.numscenarios):
                    self.buildMinBound(scenarios[s].z2, m, self.z2bounds[b], s)

            elif self.z2bounds[b] != 0 and b == 1:
                #mean constraint
                for s in range(data.numscenarios):
                    self.buildMeanBound(scenarios[s].z2, m, self.z2bounds[b], s)

            elif self.z2bounds[b] > 0 and b == 2:
                #max constraint
                for s in range(data.numscenarios):
                    self.buildMaxBound(scenarios[s].z2, m, self.z2bounds[b], s)

            elif self.z2bounds[b] > 0 and b == 3:
                for s in range(data.numscenarios):
                    self.z2eud = self.buildEUDBound(scenarios[s].z2, m, self.z2bounds[b], data, s)[0]

        print 'Generating bounds for structure number', self.index, '(', self.name, ')for zS'
        #zSbounds
        for b in range(len(self.zSbounds)):
            if self.zSbounds[b] > 0 and b == 0:
                #min constraint
                for s in range(data.numscenarios):
                    self.buildMinBound(scenarios[s].zS, m, self.zSbounds[b], s)

            elif self.zSbounds[b] != 0 and b == 1:
                #mean constraint
                for s in range(data.numscenarios):
                    self.buildMeanBound(scenarios[s].zS, m, self.zSbounds[b], s)

            elif self.zSbounds[b] > 0 and b == 2:
                #max constraint
                for s in range(data.numscenarios):
                    self.buildMaxBound(scenarios[s].zS, m, self.zSbounds[b], s)

            elif self.zSbounds[b] > 0 and b == 3:
                for s in range(data.numscenarios):
                    self.zSeud = self.buildEUDBound(scenarios[s].zS, m, self.zSbounds[b], data, s)[0]


    def buildMinBound(self, doseVector, m, bound, scen=-1):
        print "Building min bound on structure", self.index,
        if scen != -1:
            print 'scenario', scen
        else:
            print ""
        for j in range(self.size):
            doseVector[self.voxels[j]].setAttr("LB", bound)
        m.update()

    def buildMaxBound(self, doseVector, m, bound, scen=-1):
        print "Building max bound on structure", self.index,
        if scen != -1:
            print 'scenario', scen
        else:
            print ""
        for j in range(self.size):
            doseVector[self.voxels[j]].setAttr("UB", bound)
        m.update()

    def buildMeanBound(self, doseVector, m, bound, scen=-1):
        print "Building mean bound on structure", self.index,
        if scen != -1:
            print 'scenario', scen
        else:
            print ""
        if bound > 0:
            meanHolderVar = m.addVar(lb=-GRB.INFINITY, ub=bound, vtype=GRB.CONTINUOUS)
        elif bound <0:
            meanHolderVar = m.addVar(lb=-bound, vtype=GRB.CONTINUOUS)
        else:
            meanHolderVar = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
        m.update()
        m.addConstr(quicksum(doseVector[self.voxels[j]] for j in range(self.size)), GRB.EQUAL,
                    self.size * meanHolderVar, name='meanBoundConstr_' + str(self.index))
        m.update()
        return meanHolderVar

    def buildEUDBound(self, doseVector, m, bound, data, scen=-1):
        print "Building eud bound on structure", self.index,
        if scen != -1:
            print 'scenario', scen
        else:
            print ""
        # build mean holder
        meanHolderVar = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS)
        m.update()
        m.addConstr(quicksum(doseVector[self.voxels[j]] for j in range(self.size)), GRB.EQUAL, meanHolderVar)
        m.update()
        # build upper or lower bound
        boundHolderVar = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
        m.update()
        if self.index in data.targets and bound>0:
            for j in range(self.size):
                m.addConstr(boundHolderVar, GRB.LESS_EQUAL, doseVector[self.voxels[j]])
            doseEUD = m.addVar(lb=bound, vtype=GRB.CONTINUOUS)
        elif bound>0:
            for j in range(self.size):
                m.addConstr(boundHolderVar, GRB.GREATER_EQUAL, doseVector[self.voxels[j]])
            doseEUD = m.addVar(lb=0, ub=bound, vtype=GRB.CONTINUOUS)
        else:
            for j in range(self.size):
                m.addConstr(boundHolderVar, GRB.LESS_EQUAL, doseVector[self.voxels[j]])
            doseEUD = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
        m.update()
        m.addConstr(doseEUD, GRB.EQUAL, data.structGamma[self.index - 1] * meanHolderVar + (
            1 - data.structGamma[self.index - 1]) * boundHolderVar, name='eudConstr_' + str(self.index))
        m.update()
        return doseEUD, meanHolderVar



































