__author__ = 'troy'

from adaptiveClasses import *
from nonadaptiveClasses import *
from math import exp, log
from scipy import io
import numpy as np


def stochSolve(datafilename, adaptivefilename, caseloc='', alpha=-1):
    # initizlize instance of stochastic model this basically does everything other than solve the model.
    # If you want it to only add bits, then just change how the constructor works (i.e. add in
    # functions for tasks like adding the constraints, making dose, etc outside of the constructor and
    # call them from here
    print 'Starting stochastic model for alpha = ', alpha
    stoch_mod = imrt_stochastic_model(datafilename, adaptivefilename, manualAlpha=alpha, caselocation=caseloc)

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


def nonAdaSolve(datafilename, adaptivefilename, meanbound, caseloc='', ):
    nonada_mod = imrt_model_nonada(datafilename, adaptivefilename, mldbound=meanbound, caselocation=caseloc)
    nonada_mod.callSolver()
    nonada_mod.initializeCleanupNonADA()
    nonada_mod.callSolver()
    nonada_mod.outputVariablesnonADA()


def simulateAda(datafilename, adaptivefilename, dosefilenameList, biomarkerfilename, caselocation=''):
    # Goal of simulate is to calc objective functions of RILT and LTC given different biomarkers
    # For this, we load in the mean lung doses and pltc's and then evaluate these values using objective function evaluation
    # need to read in the ltc/rilt params
    data = imrt_data(datafilename, adaptivefilename, caselocation, loadDose=False)
    adaMatFile = io.loadmat(data.basedir + adaptivefilename)
    nscen = int(adaMatFile['nscen'])
    beta0 = float(adaMatFile['beta0'])
    beta1 = float(adaMatFile['beta1'])
    gamma02 = float(adaMatFile['gamma02'])
    gamma12 = float(adaMatFile['gamma12'])
    gamma22 = float(adaMatFile['gamma22'])
    s02 = float(adaMatFile['s02'])
    stageOneFraction = float(adaMatFile['s1frac'])
    biomarkers_binned = np.array(adaMatFile['biomarkers']).flatten()
    np.append(np.append([0], biomarkers_binned), [1000])
    scenprobs_binned = np.array(adaMatFile['scenprob']).flatten()

    # read in which structures go to which objective
    ptvstructindex = int(adaMatFile['ptvStruct']) - 1  # This is because of the 1 indexing in matlab
    lungstructindex = int(adaMatFile['lungStruct']) - 1
    option = int(adaMatFile['option'])

    ptvstruct = imrt_structure(data, ptvstructindex + 1)
    lungstruct = imrt_structure(data, lungstructindex + 1)



    # read in biomarkers
    bioMatFile = io.loadmat(biomarkerfilename)
    biomarkers_full = bioMatFile['fullbiomarkers']

    for dosefilename in dosefilenameList:
        # read in dose
        doseMatFile = io.loadmat(dosefilename)
        print 'opening', dosefilename
        dose1 = doseMatFile['z1']
        doseS = doseMatFile['zS']
        alpha = doseMatFile['alpha']
        doseTot = (1 - stageOneFraction) * doseS
        for i in range(doseS.shape[0]):
            doseTot[i, :] += stageOneFraction * dose1.flatten()

        Prilt = []
        Pltc = []

        for b in range(len(biomarkers_full)):
            bio = biomarkers_full[b]


            # calculate the objective given each dose

            for binIndex in range(len(biomarkers_binned) - 1):
                if bio > biomarkers_binned[binIndex] and bio <= biomarkers_binned[binIndex + 1]:
                    mld = np.mean(doseTot[binIndex, np.ix_(lungstruct.voxels)])
                    ptvmean = np.mean(doseTot[binIndex, np.ix_(ptvstruct.voxels)])
                    ptvmin = np.min(doseTot[binIndex, np.ix_(ptvstruct.voxels)])
                    # ptvmin = 0
                    ptv_eud = data.structGamma[ptvstructindex] * ptvmean + (1 - data.structGamma[
                        ptvstructindex]) * ptvmin
                    Pltc.append(getPltcFtn(ptv_eud, s02, beta0, beta1))
                    Prilt.append(getPriltFtn(mld, gamma02, gamma12, gamma22, bio))



        # output calculated Prilt and Pltc
        outFilename = data.basedir + dosefilename[:-4] + '_metrics'
        io.savemat(outFilename,
                   {'biomarkers_binned': biomarkers_binned, 'option': option, 'alpha': alpha, 'Prilt': np.array(Prilt),
                    'Pltc': np.array(Pltc), 'dosefilename': dosefilename})


def simulateNonAda(datafilename, adaptivefilename, dosefilenameList, biomarkerfilename, caselocation=''):
    # Goal of simulate is to calc objective functions of RILT and LTC given different biomarkers
    # For this, we load in the mean lung doses and pltc's and then evaluate these values using objective function evaluation
    # need to read in the ltc/rilt params
    data = imrt_data(datafilename, adaptivefilename, caselocation, loadDose=False)
    adaMatFile = io.loadmat(data.basedir + adaptivefilename)
    nscen = int(adaMatFile['nscen'])
    beta0 = float(adaMatFile['beta0'])
    beta1 = float(adaMatFile['beta1'])
    gamma02 = float(adaMatFile['gamma02'])
    gamma12 = float(adaMatFile['gamma12'])
    gamma22 = float(adaMatFile['gamma22'])
    s02 = float(adaMatFile['s02'])
    stageOneFraction = float(adaMatFile['s1frac'])
    biomarkers_binned = np.array(adaMatFile['biomarkers']).flatten()
    np.append(np.append([0], biomarkers_binned), [1000])
    scenprobs_binned = np.array(adaMatFile['scenprob']).flatten()

    # read in which structures go to which objective
    ptvstructindex = int(adaMatFile['ptvStruct']) - 1  # This is because of the 1 indexing in matlab
    lungstructindex = int(adaMatFile['lungStruct']) - 1
    option = int(adaMatFile['option'])

    ptvstruct = imrt_structure(data, ptvstructindex + 1)
    lungstruct = imrt_structure(data, lungstructindex + 1)



    # read in biomarkers
    bioMatFile = io.loadmat(biomarkerfilename)
    biomarkers_full = bioMatFile['fullbiomarkers']

    for dosefilename in dosefilenameList:
        # read in dose
        doseMatFile = io.loadmat(dosefilename)
        print 'opening', dosefilename
        dose = np.array(doseMatFile['z']).flatten()
        mldbound = doseMatFile['mldbound']
        Prilt = []
        Pltc = []

        for b in range(len(biomarkers_full)):
            bio = biomarkers_full[b]


            # calculate the objective given each dose



            mld = np.mean(dose[np.ix_(lungstruct.voxels)])
            ptvmean = np.mean(dose[np.ix_(ptvstruct.voxels)])
            ptvmin = np.min(dose[np.ix_(ptvstruct.voxels)])

            ptv_eud = data.structGamma[ptvstructindex] * ptvmean + (1 - data.structGamma[ptvstructindex]) * ptvmin
            Pltc.append(getPltcFtn(ptv_eud, s02, beta0, beta1))
            Prilt.append(getPriltFtn(mld, gamma02, gamma12, gamma22, bio))



        # output calculated Prilt and Pltc
        outFilename = data.basedir + dosefilename[:-4] + '_metrics'
        io.savemat(outFilename, {'biomarkers_binned': biomarkers_binned, 'option': option, 'mldbound': mldbound,
                                 'Prilt': np.array(Prilt), 'Pltc': np.array(Pltc), 'dosefilename': dosefilename})


def getPltcFtn(x, s02, b0, b1):
    return pow(s02, exp(b0 - b1 * x))


def getPriltFtn(x, g0, g1, g2, b):
    return (exp(g0 + (g1 + g2 * b) * x)) / (1 + exp(g0 + (g1 + g2 * b) * x))

# Function evaluations

