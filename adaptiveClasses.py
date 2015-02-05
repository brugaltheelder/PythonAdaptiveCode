__author__ = 'troy'

from scipy import io
import scipy.sparse as sps
import numpy as np
from math import exp, log
from gurobipy import *



# Data class - This contains only data and nothing else. There are no functions in the data class
class imrt_data(object):
    def __init__(self, inputFilename, adaptiveFilename):
        # This reads in the main data file
        matFile = io.loadmat(inputFilename)
        #number voxels, number bixels, num nonzeros in D matrix, number structs
        self.nVox, self.nBix, self.nDijs, self.nStructs = int(matFile['nvox']), int(matFile['nbixel']), int(
            matFile['numdijs']), int(matFile['numstructs'])
        #number oas, number targets, list of oar indices, list of target indices
        self.numoars, self.numtargets, self.oars, self.targets = int(matFile['numoars']), int(
            matFile['numtargets']), np.array(matFile['oars']).flatten(), np.array(matFile['targets']).flatten()

        # These are holders for bixel and voxel indices and dijs values. They assume indexing starts at 1 (matlab)
        #NOTE: You need to run my conversion script in MATLAB to get the proper dose delivered (rebuildVoxelIndices.m)
        bixe = np.array(matFile['bixe2_new']).flatten() - 1
        voxe = np.array(matFile['voxe2_new_nvox']).flatten() - 1
        dijs = np.array(matFile['dijs2_new']).flatten()

        # Build the sparse matrix
        self.Dmat = sps.csr_matrix((dijs, (bixe, voxe)), shape=(self.nBix, self.nVox))
        #Reads in mask value (unused due to "structs" matlab array
        self.maskValue = np.array(matFile['maskValue']).flatten()
        #Structure index (starting at 1) per voxel
        self.structPerVoxel = np.array(matFile['structs']).flatten()
        # Names of organs, mapped to strings
        self.pickstructs = map(str, matFile['pickstructs'])



        # Read in structure bounds matrix: format is:
        # Each row is a structure
        # Each column is the following:
        # col0-3: z1minbound, z1meanbound, z1maxbound, z1eudbound
        # col4-7: z2minbound, z2meanbound, z2maxbound, z2eudbound
        # col8-11: zSminbound, zSmeanbound, zSmaxbound, zSeudbound
        self.structBounds = np.array(matFile['structurebounds'])
        # This is the EUD weight for each structure
        self.structGamma = np.array(matFile['eudweights']).flatten()

        # Read in the filename for your particular method-specifi class
        self.adaptiveFilename = adaptiveFilename
        # Open the file
        adaptiveFile = io.loadmat(self.adaptiveFilename)
        # Determine how many scenarios here. This variable has to be in your matlab file. Same thing with s1frac
        self.numscenarios = int(adaptiveFile['nscen'])
        self.stageOneFraction = float(adaptiveFile['s1frac'])


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
            self.structures[s].buildConstraints(self.data, self.m, self.z1, self.scenarios)

        # THIS IS THE IMPORTANT PART - add in your class here and remove mine. Mine, upon construction, builds the necessary variables and constraints
        # necessary for my method. I'm going to add a few other inputs to it so I can mass-run things, but those will come later
        # Use this as a template
        print 'initializing adaptive lung class'
        self.adaLung = imrt_adaptiveLung(self.data, self.m, self.scenarios, self.structures, self.z1)

        # Uncomment to write out model
        # print 'Writing out model'
        # self.m.write('out.lp')
        #print 'Model writing done'

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


# Structure class - This holds the structure-specific bounds and associated constraints
class imrt_structure(object):
    def __init__(self, data, index):
        assert (isinstance(data, imrt_data))
        self.name = data.pickstructs[
            index - 1]  # ASSUMES ANAT INDEXING STARTS AT 1 TODO FIX THIS SO IT STARTS AT 0, also below and elsewhere
        self.index = index
        self.voxels = np.where(data.structPerVoxel == index)[0]
        self.size = self.voxels.size
        print self.index, self.size, self.name
        # This gets the bounds from structBounds (see imrt_data comment on structBounds for details on structure)
        self.z1bounds = data.structBounds[index - 1, 0:4]
        self.z2bounds = data.structBounds[index - 1, 4:8]
        self.zSbounds = data.structBounds[index - 1, 8:12]

    # This goes through each given bound for z1, z2, zS for this structure and builds the constraint as necessary
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


    # Sets lower bound on dose to each voxel in the structure
    def buildMinBound(self, doseVector, m, bound, scen=-1):
        print "Building min bound on structure", self.index,
        if scen != -1:
            print 'scenario', scen
        else:
            print ""
        for j in range(self.size):
            doseVector[self.voxels[j]].setAttr("LB", bound)
        m.update()

    # Sets upper bound on dose to each voxel in the structure
    def buildMaxBound(self, doseVector, m, bound, scen=-1):
        print "Building max bound on structure", self.index,
        if scen != -1:
            print 'scenario', scen
        else:
            print ""
        for j in range(self.size):
            doseVector[self.voxels[j]].setAttr("UB", bound)
        m.update()

    # Builds a mean dose variable and returns that gurobi variable (just in case you want it)
    # Ignore the scen = -1...I wanted to do something cool with the naming, but gave up.
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
            #This is included just in case you don't want a bounded mean variable
            meanHolderVar = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
        m.update()
        m.addConstr(quicksum(doseVector[self.voxels[j]] for j in range(self.size)), GRB.EQUAL,
                    self.size * meanHolderVar, name='meanBoundConstr_' + str(self.index) + self.name.strip())
        m.update()
        return meanHolderVar


    # This function builds and EUD. It assumes you want max and mean for OAR and min and mean for Target
    # This will also generate unbouned EUD vars (set input bound to zero). I also have a sanity check "makeIfZero",
    # which NEEDS TO BE TRUE IF YOU WANT TO GENERATE AN UNBOUNDED EUD
    def buildEUDBound(self, doseVector, m, bound, data, scen=-1, makeIfZero=False):
        print "Building eud bound on structure", self.index,
        if scen != -1:
            print 'scenario', scen
        else:
            print ""
        # build mean holder
        meanHolderVar = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS)
        m.update()
        m.addConstr(quicksum(doseVector[self.voxels[j]] for j in range(self.size)), GRB.EQUAL, self.size * meanHolderVar)
        m.update()
        # build upper or lower bound depending on oar or target
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
        elif makeIfZero and self.index in data.targets:
            for j in range(self.size):
                m.addConstr(boundHolderVar, GRB.LESS_EQUAL, doseVector[self.voxels[j]])
            doseEUD = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
        elif makeIfZero and self.index in data.oars:
            for j in range(self.size):
                m.addConstr(boundHolderVar, GRB.GREATER_EQUAL, doseVector[self.voxels[j]])
            doseEUD = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
        m.update()
        m.addConstr(doseEUD, GRB.EQUAL, data.structGamma[self.index - 1] * meanHolderVar + (
            1 - data.structGamma[self.index - 1]) * boundHolderVar,
                    name='eudConstr_' + str(self.index) + self.name.strip())
        m.update()
        return doseEUD, meanHolderVar



































