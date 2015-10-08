'''
Non-profiled LRA a la ASIACRYPT'13 paper [https://eprint.iacr.org/2013/794].
The target is AES encryption with the S-box output as the target variable.
Optional pre-averaging of traces and normalization for goodness-of-fit
traces are implemented. CPA is included for comparison.

Uses manual OLS (dot-products and matrix inversion, relying on numpy-MKL
efficient implementation).

This file is a library of functions. The attacks are supposed to be
implemented as scripts importing this library.

Heavily under development, correctness not guaranteed.

Version: 0.1, 2014-11-19
Started by Ilya on 2014-05-25
'''

import numpy as np

# preload the precomputed lookup tables, to avoid bloating of this code
sbox              = np.load('data/aessbox.npy')           # AES S-box
invsbox           = np.load('data/aesinvsbox.npy')        # AES inverse S-box
byteHammingWeight = np.load('data/bytehammingweight.npy') # HW of a byte


#############################################################################
### Common for both attacks

# averaging of traces based on the value of the data byte
# TODO: remove, implemented in condaver.py
def conditionalAveragingAESSbox(data, traces):

    (numTraces, traceLength) = traces.shape

    # array for averaged traces
    avtraces = np.zeros((256, traceLength)) # TODO check if need finer precision
    counters = np.zeros(256)

    # incremental averaging
    for i in xrange(numTraces):
        if (counters[data[i]]==0):
            avtraces[data[i]] = traces[i]
        else:
            avtraces[data[i]] = avtraces[data[i]] + (traces[i] - avtraces[data[i]]) / counters[data[i]]
        counters[data[i]] += 1

    # get an vector of only _observed_ observed values
    avdata = np.flatnonzero(counters)

    # remove lines corresponding to non-observed values
    avtraces = avtraces[avdata]

    return avdata, avtraces

# Functions for computing intermediate values
# data is a 1-D array, keyByte is a scalar
def sBoxOut(data, keyByte):
    sBoxIn = data ^ keyByte
    return sbox[sBoxIn]
def sBoxInXorOut(data, keyByte):
    sBoxIn = data ^ keyByte
    return sBoxIn ^ sbox[sBoxIn]
def invSboxOut(data, keyByte):
    sBoxIn = data ^ keyByte
    return invsbox[sBoxIn]
def invSboxInXorOut(data, keyByte):
    sBoxIn = data ^ keyByte
    return sBoxIn ^ invsbox[sBoxIn]

##############################################################################
### A. LRA attack stuff

### Leakge modelling
# These functions do 2 things in the same place:
# 1. define basis functions gi(x) of a leakage model for a byte x:
#    b0 x g0(x) + b1 x g1(x) + ... + bn x gn(x)
# 2. compute and return the values of gi(x), such that they can be used
#    later to obtain rows of the matrix for linear regression
# Note that column of ones is included!
# Note also that not all the functions are currently compatible with the code
#  in the CPA and LRA functions because latter use wrapper functions for
#  incremental binding. Only those are compatible that take a second bitWidth
#  argument.
# TODO: make all functions compatible with incremental binding.


# A simple linear model - sum of bits with different coefficients:
#  gi = xi, 0 <= i < bitWidth.
def basisModelSingleBits(x, bitWidth):
    g = []
    for i in range(0, bitWidth):
        bit = (x >> i) & 1  # this is the definition: gi = [bit i of x]
        g.append(bit)
    g.append(1)
    return g

# Invididual bits and all pairwise products of bits
def basisModelSingleBitsAndPairs(x, bitWidth):
    g = []
    for i in range(0, bitWidth):
        # append single bits
        bit = (x >> i) & 1  # this is the definition: gi = [bit i of x]
        g.append(bit)
        # append pairs
        for j in range(i + 1, bitWidth):
            otherbit = (x >> j) & 1
            bitproduct = bit * otherbit
            g.append(bitproduct)
    g.append(1)
    return g

# A Hamming weight model: g0 = HW(x)
def basisModelHW(x):
    g = []
    hw = byteHammingWeight[x]  # this is the definition: gi = HW(x)
    g.append(hw)
    g.append(1)
    return g

# An 'all 256 bit combinations' model:
# a) helper from http://wiki.python.org/moin/BitManipulation
def parityOf(int_type):
    parity = 0
    while (int_type):
        parity = ~parity
        int_type = int_type & (int_type - 1)
    if (parity != 0): # to convert -1 to 1
        parity = 1
    return(parity)
# b) the model itself
def basisModel256(x):
    g = []
    # note that we start from 1 to exclude case 0 which means the function
    # does not depend on any bit of x, i.e. a constant - we will add the
    # constant explicitly later as the last column.
    for i in np.arange(1, 256, dtype='uint8'):
        xmasked = x & i
        gi = parityOf(xmasked)
        g.append(gi)
    g.append(1)
    return g

# LRA attack on AES
# data                 - 1-D array of input bytes
# traces               - 2-D array of traces
# intermediateFunction - one of the functions like sBoxOut above in the common section 
# basisFunctionsModel  - one of the functions like basisModelSingleBits above
#                        in this section
# TODO parametrize hard-coded values such as 256, 8, refactor to merge common part with DES
def lraAES(data, traces, intermediateFunction, basisFunctionsModel):

    ### 0. some helper variables
    (numTraces, traceLength) = traces.shape
    
    # define a wrapper for currying (incremental parameter binding)
    def basisFunctionsModelWrapper(y):
        def basisFunctionsModelCurry(x):
            return basisFunctionsModel(x, y)
        return basisFunctionsModelCurry

    ### 1: compute SST over the traces
    SStot = np.sum((traces - np.mean(traces, 0)) ** 2, 0)

    ### 2. The main attack loop

    # preallocate arrays
    SSreg = np.empty((256, traceLength)) # Sum of Squares due to regression
    E = np.empty(numTraces)              # expected values

    allCoefs = [] # placeholder for regression coefficients

    # per-keycandidate loop
    for k in np.arange(0, 256, dtype='uint8'):

        # predict intermediate variable
        intermediateVariable = intermediateFunction(data, k)

        # buld equation system
        M = np.array(map(basisFunctionsModelWrapper(8), intermediateVariable))

        # some precomputations before the per-sample loop
        P = np.dot(np.linalg.inv(np.dot(M.T, M)), M.T)
        #Q = np.dot(M, P)

        coefs = [] # placeholder for regression coefficients

        # per-sample loop: solve the system for each time moment
        for u in range(0,traceLength):

            # if do not need coefficients beta - use precomputed value
            #np.dot(Q, traces[:,u], out=E)

            # if need the coefficients - do the multiplication using
            # two dot products and let the functuion return beta alongside R2
            beta = np.dot(P, traces[:,u])
            coefs.append(beta)
            E = np.dot(M, beta)

            SSreg[k,u] = np.sum((E - traces[:,u]) ** 2)

        allCoefs.append(coefs)
        #print 'Done with candidate', k

    ### 3. compute Rsquared
    R2 = 1 - SSreg / SStot[None, :]

    return R2, allCoefs

# LRA attack on DES
# data                 - array of inputs (format depends on intermediateFunction)
# traces               - 2-D array of traces
# intermediateFunction - one of functions like sBoxOut above in the common section
# sBoxNumber           - DES S-box to attack
# basisFunctionsModel  - one of function like basisModel9 above in this section
# TODO parametrize hard-coded values such as 64, 6, refactor to merge common part with AES
def lraDES(data, traces, intermediateFunction, sBoxNumber, basisFunctionsModel):

    ### 0. some helper variables
    (numTraces, traceLength) = traces.shape

    # define a wrapper (currying) for incremental parameter binding
    def basisFunctionsModelWrapper(y):
        def basisFunctionsModelCurry(x):
            return basisFunctionsModel(x, y)
        return basisFunctionsModelCurry

    ### 1: compute SST over the traces
    SStot = np.sum((traces - np.mean(traces, 0)) ** 2, 0)

    ### 2. The main attack loop

    # preallocate arrays
    SSreg = np.empty((64, traceLength)) # Sum of Squares due to regression
    E = np.empty(numTraces)              # expected values

    # per-keycandidate loop
    for k in np.arange(0, 64, dtype='uint8'):

        # predict intermediate variable
        intermediateVariable = intermediateFunction(data, k, sBoxNumber)

        # buld equation system
        M = np.array(map(basisFunctionsModelWrapper(6), intermediateVariable))

        # some precomputations before the per-sample loop
        P = np.dot(np.linalg.inv(np.dot(M.T, M)), M.T)
        Q = np.dot(M, P)

        # per-sample loop: solve the system for each time moment
        for u in range(0,traceLength):

            # if do not need coefficients beta - use precomputed value
            np.dot(Q, traces[:,u], out=E)

            # if need the coefficients - do the multiplication using
            # two dot products and let the functuion return beta alongside R2
            #beta = np.dot(P, traces[:,u])
            #E = np.dot(M, beta)

            SSreg[k,u] = np.sum((E - traces[:,u]) ** 2)

        #print 'Done with candidate', k

    ### 3. compute Rsquared
    R2 = 1 - SSreg / SStot[None, :]

    return R2

# convert R2 to adjusted R2 (https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2)
# n - number of samples
# p - the total number of regressors in the linear model (i.e. basis functions), excluding the linear term
def adjustedR2(R2, n, p):
    R2adj = 1 - ((1 - R2 ** 2) * (n - 1) / np.double(n - p - 1))
    return R2adj

# normalize the matrix of distinguisher traces according to ASIACRYPT'13 proposal
def normalizeR2Traces(R2):
    R2norm = np.empty(R2.shape)
    traceLength = R2.shape[1]
    for i in range(0,traceLength): # TODO should be possible to do it in one line without a loop
        R2norm[:,i] = (R2[:,i] - np.mean(R2[:,i])) / np.var(R2[:,i])
    return R2norm

##############################################################################
### A. CPA attack stuff

def leakageModelHW(x):
    return byteHammingWeight[x]

# correlation trace computation as improved by StackOverflow community
# O - matrix of observed leakage (i.e. traces)
# P - column of predictions
# returns a correlation trace
def correlationTraceSO(O, P):
    n = P.size
    DO = O - (np.einsum('ij->j', O, dtype='float64') / np.double(n))
    DP = P - (np.einsum('i->', P, dtype='float64') / np.double(n))
    tmp = np.einsum('ij,ij->j', DO, DO)
    tmp *= np.einsum('i,i->', DP, DP)
    return np.dot(DP, DO) / np.sqrt(tmp)

# CPA attack
# data                 - 1-D array of input bytes
# traces               - 2-D array of traces
# intermediateFunction - one of functions like sBoxOut above in the common section
# leakageFunction      - one of the fucntions like leakgeModelHW above in this section
def cpaAES(data, traces, intermediateFunction, leakageFunction):

    traceLength = traces.shape[1]

    # compute intermediate variable predictions
    k = np.arange(0,256, dtype='uint8') # key chunk candidates
    H = np.zeros((256, len(data)), dtype='uint8') # intermediate variable predictions
    for i in range(256):
        H[i,:] = intermediateFunction(data, k[i])

    # compute leakage hypotheses for every  all the key candidates
    HL = map(leakageFunction, H) # leakage model here (HW for now)

    CorrTraces = np.empty([256, traceLength]);

    # per-keycandidate loop
    # TODO: loop can be avoided by clever use of np.einsum in correlationTraceSO(...)
    for i in range(0, 256):
        CorrTraces[i] = correlationTraceSO(traces, HL[i])

    return CorrTraces

# CPA attack on DES
# data                 - array of inputs (format depends on intermediateFunction)
# traces               - 2-D array of traces
# intermediateFunction - one of functions like sBoxOut above in the common section
# sBoxNumber           - DES S-box to attack
# leakageFunction      - one of the fucntions like leakgeModelHW above in this section
def cpaDESwithAveraging(data, traces, intermediateFunction, sBoxNumber, leakageFunction):

    traceLength = traces.shape[1]

    # compute intermediate variable predictions
    k = np.arange(0,64, dtype='uint8') # key chunk candidates
    H = np.zeros((64, len(data)), dtype='uint8') # intermediate variable predictions
    for i in range(64):
        for j in range(0, len(data)):
            H[i,j] = intermediateFunction(data[j], k[i], sBoxNumber)

    # compute leakage hypotheses for every  all the key candidates
    HL = map(leakageFunction, H) # leakage model here (HW for now)

    CorrTraces = np.empty([64, traceLength]);

    # per-keycandidate loop
    # TODO: loop can be avoided by clever use of np.einsum in correlationTraceSO(...)
    for i in range(0, 64):
        CorrTraces[i] = correlationTraceSO(traces, HL[i])

    return CorrTraces
