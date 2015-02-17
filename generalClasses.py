__author__ = 'troy'
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