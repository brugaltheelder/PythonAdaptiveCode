__author__ = 'troy'

from scipy import io
import scipy.sparse as sps
import numpy as np
from gurobipy import *



# Data class
class imrt_data(object):
    def __init__(self, inputFilename):
        matFile = io.loadmat(inputFilename)
        self.nVox, self.nBix, self.nDijs, self.nStructs = int(matFile['nvox']), int(matFile['nbixel']), int(
            matFile['numdijs']), int(matFile['numstructs'])
        self.numoars, self.numtargets, self.oars, self.targets = int(matFile['numoars']), int(
            matFile['numtargets']), np.array(matFile['oars']).flatten(), np.array(matFile['targets']).flatten()
        print self.nVox, self.nBix, self.nDijs, self.oars
        bixe = np.array(matFile['bixe2_new']).flatten() - 1
        voxe = np.array(matFile['voxe2_new_nvox']).flatten() - 1
        dijs = np.array(matFile['dijs2_new']).flatten()
        self.Dmat = sps.csr_matrix((dijs, (bixe, voxe)), shape=(self.nBix, self.nVox))
        self.maskValue = np.array(matFile['maskValue']).flatten()
        self.structPerVoxel = np.array(matFile['structs']).flatten()
        self.pickstructs = map(str, matFile['pickstructs'])

        self.structBounds = np.array(matFile['structurebounds'])
        self.structGamma = np.array(matFile['eudweights']).flatten()

        # todo read in scenario data file location and data file
        self.numscenarios = 2
        self.scneariovalues = numpy.array([1.2, 2.3])


class scenario(object):
    def __init__(self, data, num, m):
        assert (isinstance(data, imrt_data))
        self.num = num
        self.scenValue = data.scneariovalues(self.num)
        print 'building scenario', self.num
        # build dose variables
        self.z2 = [m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in xrange(data.nVox)]
        m.update()
        #initialize dose constraint
        self.doseConstr2 = [m.addConstr(-self.z2[j], GRB.EQUAL, 0) for j in xrange(data.nVox)]
        m.update()
        #add in beamlet intensities
        self.x2 = [m.addVar(lb=0, vtype=GRB.CONTINUOUS,
                            column=Column(np.array(data.Dmat.getrow(i).todense()).flatten().tolist(), self.doseConstr2))
                   for i in xrange(data.nBix)]
        m.update()

        # todo




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
        self.z1 = [self.m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in xrange(self.data.nVox)]
        self.m.update()
        print 'Building Stage-one Dose Constraint'
        self.doseConstr1 = [self.m.addConstr(-self.z1[j], GRB.EQUAL, 0) for j in xrange(self.data.nVox)]
        self.m.update()
        print 'Populating Stage-one Dose Constraint'
        self.x1 = [self.m.addVar(lb=0, vtype=GRB.CONTINUOUS,
                                 column=Column(np.array(self.data.Dmat.getrow(i).todense()).flatten().tolist(),
                                               self.doseConstr1)) for i in xrange(self.data.nBix)]
        self.m.update()
        print 'Stage-one Dose Constraint Built'



        # Uncomment to write out model
        # print 'Writing out model'
        #self.m.write('out.lp')
        #print 'Model writing done'

        #todo Initialize scenarios (which build other gurobi variables for overall Zs)
        self.scenarios = [scenario(self.data, s, self.m) for s in range(self.data.numscenarios)]

        #todo initizlize stochastic class

    def callSolver(self):
        self.m.optimize()

        # todo make function that builds bounds for each structure (sets of bounds: Z1, Z2S, ZS, min mean max eud)


# Structure class
class imrt_structure(object):
    def __init__(self, data, index):
        assert (isinstance(data, imrt_data))
        self.name = data.pickstructs[
            index - 1]  # ASSUMES ANAT INDEXING STARTS AT 1 TODO FIX THIS SO IT STARTS AT 0, also below
        self.index = index
        self.voxels = np.where(data.structPerVoxel == index)[0]
        self.size = self.voxels.size
        self.z1bounds = data.structBounds[index - 1, 0:4]
        self.z2bounds = data.structBounds[index - 1, 4:8]
        self.zSbounds = data.structBounds[index - 1, 8:12]



# adaLung class

# todo reads in separate data file


