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
        # number voxels, number bixels, num nonzeros in D matrix, number structs
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
        #Reads in mask value (unused due to "structs" matlab array)
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
        # This is the EUD weight (for mean) for each structure
        self.structGamma = np.array(matFile['eudweights']).flatten()

        # Read in the filename for your particular method-specifi class
        self.adaptiveFilename = adaptiveFilename
        # Open the file
        adaptiveFile = io.loadmat(self.adaptiveFilename)
        # Determine how many scenarios here. This variable has to be in your matlab file. Same thing with s1frac
        self.numscenarios = int(adaptiveFile['nscen'])
        self.stageOneFraction = float(adaptiveFile['s1frac'])


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


    def buildConstraintsNonAdaptive(self, data, m, zdose):
        # for each set of bounds, for each bound value (if >0), build constraint
        print 'Generating bounds for structure number', self.index, '(', self.name, ')for z'
        # z1bounds
        for b in range(len(self.z1bounds)):
            if self.zSbounds[b] > 0 and b == 0:
                # min constraint
                self.buildMinBound(zdose, m, self.zSbounds[b])

            elif self.zSbounds[b] != 0 and b == 1:
                # mean constraint
                self.buildMeanBound(zdose, m, self.zSbounds[b])

            elif self.zSbounds[b] > 0 and b == 2:
                # max constraint
                self.buildMaxBound(zdose, m, self.zSbounds[b])

            elif self.zSbounds[b] > 0 and b == 3:
                self.zeud = self.buildEUDBound(zdose, m, self.z1bounds[b], data)[0]


    # This goes through each given bound for z1, z2, zS for this structure and builds the constraint as necessary
    def buildConstraintsAdaptive(self, data, m, z1dose, scenarios):
        # for each set of bounds, for each bound value (if >0), build constraint
        print 'Generating bounds for structure number', self.index, '(', self.name, ')for z1'
        # z1bounds
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
                self.z2eud = []
                for s in range(data.numscenarios):
                    self.z2eud.append(self.buildEUDBound(scenarios[s].z2, m, self.z2bounds[b], data, s)[0])

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
        elif bound < 0:
            meanHolderVar = m.addVar(lb=-bound, vtype=GRB.CONTINUOUS)
        else:
            # This is included just in case you don't want a bounded mean variable
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
        m.addConstr(quicksum(doseVector[self.voxels[j]] for j in range(self.size)), GRB.EQUAL,
                    self.size * meanHolderVar)
        m.update()
        # build upper or lower bound depending on oar or target
        boundHolderVar = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
        m.update()
        if self.index in data.targets and bound > 0:
            for j in range(self.size):
                m.addConstr(boundHolderVar, GRB.LESS_EQUAL, doseVector[self.voxels[j]])
            doseEUD = m.addVar(lb=bound, vtype=GRB.CONTINUOUS)
        elif bound > 0:
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








