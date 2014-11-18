'''
Non-profiled LRA a la ASIACRYPT'13 paper [https://eprint.iacr.org/2013/794].
The target is AES encryption with the S-box output as the target variable.
Optional pre-averaging of traces and normalization for goodness-of-fit
traces are implemented. CPA is included for comparison.

Uses manual OLS (dot-products and matrix inversion, relying on numpy-MKL
efficient implementation).

The following results were obtained on an Intel Core i5-2540M, 8GB RAM,
Windows 7 x64, Python 2.7.8 x64, numpy-MKL 1.9.0.

In [3]: %run lra
Traceset parameters
Number of traces: 2000
Trace length: 200
---
Attacks with 1000 traces:
Running CPA... done in 0.349224 s
Running LRA... done in 37.254000 s
Normalizing LRA results... done
---
Attacks with 1000 traces and conditional averaging:
Performing conditional trace averaging... done in 0.003820 s
Running CPA on averaged traces... done in 0.057802 s
Running LRA on averaged traces... done in 3.997455 s
Normalizing LRA results... done
---
Plotting...

Started by Ilya on 2014-05-25
'''

import numpy as np
import matplotlib.pyplot as plt
import time

#############################################################################
### Trace averaging: common for both attacks

# averaging of traces based on the value of the data byte
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

##############################################################################
### A. LRA attack stuff

### Leakge modelling
# These functions do 2 things in the same place:
# 1. define basis functions gi(x) of a leakage model for a byte x:
#    b0 x g0(x) + b1 x g1(x) + ... + bn x gn(x)
# 2. compute and return the values of gi(x), such that they can be used
#    later to obtain rows of the matrix for linear regression
# Note that column of ones is included!

# A simple 9-component linear model (sum of bits with different
#  coefficients): gi = xi, 0 <= i < 8.
def leakageModel9(x):
    g = []
    for i in range(0, 8):
        bit = (x >> i) & 1  # this is the definition: gi = [bit i of x]
        g.append(bit)
    g.append(1)
    return g
# same but only two LSB's included (because teh above shows that they are the
#  most contributing ones)
def leakageModel2(x):
    g = []
    for i in range(0, 2):
        bit = (x >> i) & 1
        g.append(bit)
    g.append(1)
    return g

# A Hamming weight model: g0 = HW(x)
def leakageModelHW(x):
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
def leakageModel256(x):
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
def lraAES(data, traces):

    (numTraces, traceLength) = traces.shape

    ### 1: compute SST over the traces
    SStot = np.sum((traces - np.mean(traces, 0)) ** 2, 0)

    ### 2. The main attack loop

    # preallocate arrays
    SSreg = np.empty((256, traceLength)) # Sum of Squares due to regression
    E = np.empty(numTraces)              # expected values

    # per-keycandidate loop
    for k in np.arange(0, 256, dtype='uint8'):

        # predict intermediate variable
        sBoxOut = sbox[data ^ k]

        # buld equation system
        M = np.array(map(leakageModel9, sBoxOut))

        # some precomputations before the per-sample loop
        P = np.dot(np.linalg.inv(np.dot(M.T, M)), M.T)
        Q = np.dot(M, P)

        # per-sample loop: solve the system for each time moment
        for u in range(0,traceLength):

            # if do not need coefficients beta - use precomputed value
            np.dot(Q, traces[:,u], out=E)

            # if need the coefficients - do the multiplication using
            # two dot products and let the functuion return beta alongside R2
            #beta = np.dot(P, avtraces_observed[:,u])
            #E = np.dot(M, beta)

            SSreg[k,u] = np.sum((E - traces[:,u]) ** 2)

        #print 'Done with candidate', k

    ### 3. compute Rsquared
    R2 = 1 - SSreg / SStot[None, :]

    return R2

# normalize the matrix of distinguisher traces according to ASIACRYPT'13 proposal
def normalizeR2Traces(R2):
    R2norm = np.empty(R2.shape)
    traceLength = R2.shape[1]
    for i in range(0,traceLength): # TODO should be possible to do it in one line without a loop
        R2norm[:,i] = (R2[:,i] - np.mean(R2[:,i])) / np.var(R2[:,i])
    return R2norm

##############################################################################
### A. CPA attack stuff

# preload the precomputed lookup tables
byteHammingWeight = np.load('data/bytehammingweight.npy') # HW of a byte
sbox = np.load('data/aessbox.npy') # AES S-box

def leakageModel_HW(x):
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
def cpaAES(data, traces):

    traceLength = traces.shape[1]

    k = np.arange(0,256, dtype='uint8') # key chunk candidates
    H = np.zeros((256, len(data)), dtype='uint8') # intermediate variable predictions
    for i in range(256):
        H[i,:] = sbox[data ^ k[i]]

    # compute leakage hypotheses for every  all the key candidates
    HL = map(leakageModel_HW, H) # leakage model here (HW for now)

    CorrTraces = np.empty([256, traceLength]);

    # per-keycandidate loop
    for i in range(0, 256):
        CorrTraces[i] = correlationTraceSO(traces, HL[i])

    return CorrTraces


##############################################################################
### Different test flows (TODO: move all below this line to another file,
###  leaving this file a library only, with a couple of tests)

def testComparePerformance():

    ### 0. Prerequisites

    # the correct key byte (for later metrics)
    correctKey = np.array([0x2B,0x7E,0x15,0x16,0x28,0xAE,0xD2,0xA6,0xAB,0xF7,0x15,0x88,0x09,0xCF,0x4F,0x3C])

    ### 1. Load samples and data

    # parameters
    SboxNum = 0

    # readout
    npzfile = np.load('traces/swaes_atmega_powertraces.npz')
    data = npzfile['data'][:,SboxNum]      # selecting only the required byte
    traces = npzfile['traces'][:,900:1100] # selecting the sub-range of samples

    # output traceset parameters
    (numTraces, traceLength) = traces.shape
    print "Traceset parameters"
    print "Number of traces:", numTraces
    print "Trace length:", traceLength

    ### 2. LRA and CPA with a fixed number of traces

    N = 200 # number of traces

    print "---\nAttacks with %d traces" % N
    print "Running CPA...",
    t0 = time.clock()
    CorrTraces = cpaAES(data[0:N], traces[0:N])
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Running LRA...",
    t0 = time.clock()
    R2 = lraAES(data[0:N], traces[0:N])
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Normalizing LRA results...",
    R2norm = normalizeR2Traces(R2)
    print "done"

    print "---\nAttacks with %d traces and conditional averaging" % N
    print "Performing conditional trace averaging...",
    t0 = time.clock()
    (avdata, avtraces) = conditionalAveragingAESSbox(data[0:N], traces[0:N])
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Running CPA on averaged traces...",
    t0 = time.clock()
    CorrTracesAv = cpaAES(avdata, avtraces)
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Running LRA on averaged traces...",
    t0 = time.clock()
    R2Av = lraAES(avdata, avtraces)
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Normalizing LRA results...",
    R2AvNorm = normalizeR2Traces(R2Av)
    print "done"

    ### 3.  visualize the result, highlighting the correct trace

    print "---\nPlotting..."
    fig, ax = plt.subplots(3,2,sharex=True, squeeze=True)

    WrongKeyRange = range(0, correctKey[SboxNum]) + range(correctKey[SboxNum] + 1, 256)

    ax[0][0].plot(CorrTraces[WrongKeyRange, :].T, color = 'grey')
    ax[0][0].plot(CorrTraces[correctKey[SboxNum], :], 'r')

    ax[1][0].plot(R2[WrongKeyRange, :].T, color = 'grey')
    ax[1][0].plot(R2[correctKey[SboxNum], :], 'r')

    ax[2][0].plot(R2norm[WrongKeyRange, :].T, color = 'grey')
    ax[2][0].plot(R2norm[correctKey[SboxNum], :], 'r')

    ax[0][1].plot(CorrTracesAv[WrongKeyRange, :].T, color = 'grey')
    ax[0][1].plot(CorrTracesAv[correctKey[SboxNum], :], 'r')

    ax[1][1].plot(R2Av[WrongKeyRange, :].T, color = 'grey')
    ax[1][1].plot(R2Av[correctKey[SboxNum], :], 'r')

    ax[2][1].plot(R2AvNorm[WrongKeyRange, :].T, color = 'grey')
    ax[2][1].plot(R2AvNorm[correctKey[SboxNum], :], 'r')

    # same vertical scales for correlation and R2
    ax[0][0].set_ylim(ax[0][1].get_ylim())
    ax[1][0].set_ylim(ax[1][1].get_ylim())

    fig.suptitle("CPA and LRA on %d traces" % N)
    ax[0][0].set_title('Without cond. averaging')
    ax[0][1].set_title('With cond. averaging')
    ax[0][0].set_ylabel('Correlation')
    ax[1][0].set_ylabel('R2')
    ax[2][0].set_ylabel('Normalized R2')
    ax[2][0].set_xlabel('Time sample')
    ax[2][1].set_xlabel('Time sample')

    plt.show()

def testSeriosAttack():

    ### 0. Prerequisites

    # the correct key byte (for later metrics)
    correctKey = np.array([0x2B,0x7E,0x15,0x16,0x28,0xAE,0xD2,0xA6,0xAB,0xF7,0x15,0x88,0x09,0xCF,0x4F,0x3C])

    ### 1. Load samples and data

    # parameters
    SboxNum = 0

    # readout
    npzfile = np.load('traces/swaes_atmega_powertraces.npz')
    data = npzfile['data'][:,SboxNum]      # selecting only the required byte
    traces = npzfile['traces'][:,900:1100] # selecting the sub-range of samples

    # output traceset parameters
    (numTraces, traceLength) = traces.shape
    print "Traceset parameters"
    print "Number of traces:", numTraces
    print "Trace length:", traceLength

    ### 2. LRA and CPA with a fixed number of traces

    N = 200 # number of traces

    print "---\nAttacks with %d traces and conditional averaging" % N
    print "Performing conditional trace averaging...",
    t0 = time.clock()
    (avdata, avtraces) = conditionalAveragingAESSbox(data[0:N], traces[0:N])
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Running CPA on averaged traces...",
    t0 = time.clock()
    CorrTracesAv = cpaAES(avdata, avtraces)
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Running LRA on averaged traces...",
    t0 = time.clock()
    R2Av = lraAES(avdata, avtraces)
    t1 = time.clock()
    print "done in %f s" % (t1 - t0)
    print "Normalizing LRA results...",
    R2AvNorm = normalizeR2Traces(R2Av)
    print "done"

    ### 3.  visualize the result, highlighting the correct trace

    print "---\nPlotting..."
    fig, ax = plt.subplots(3,1,sharex=True, squeeze=True)

    WrongKeyRange = range(0, correctKey[SboxNum]) + range(correctKey[SboxNum] + 1, 256)

    ax[0].plot(CorrTracesAv[WrongKeyRange, :].T, color = 'grey')
    ax[0].plot(CorrTracesAv[correctKey[SboxNum], :], 'r')

    ax[1].plot(R2Av[WrongKeyRange, :].T, color = 'grey')
    ax[1].plot(R2Av[correctKey[SboxNum], :], 'r')

    ax[2].plot(R2AvNorm[WrongKeyRange, :].T, color = 'grey')
    ax[2].plot(R2AvNorm[correctKey[SboxNum], :], 'r')

    # same vertical scales for correlation and R2
    ax[0].set_ylim(ax[0].get_ylim())
    ax[1].set_ylim(ax[1].get_ylim())

    fig.suptitle("CPA and LRA on %d traces" % N)
    ax[0].set_title('With cond. averaging')
    ax[0].set_ylabel('Correlation')
    ax[1].set_ylabel('R2')
    ax[2].set_ylabel('Normalized R2')
    ax[2].set_xlabel('Time sample')

    plt.show()

 
##############################################################################
### The main execution body

if __name__ == '__main__':

    # TODO: common file-reading code should be here

    # run different scenarious defined as functions above

    #testComparePerformance()

    testSeriosAttack()


