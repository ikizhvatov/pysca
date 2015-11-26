'''
Compare perfromance of LRA and CPA with and without conditional averaging

pysca toolbox
Version: 0.3, 2015-10-22
Started by Ilya on 2014-11-18
'''

import numpy as np
import matplotlib.pyplot as plt
import time

from aes import AES    # interweb's SlowAES toolbox
from lracpa import *   # my LRA-CPA toolbox
from condaveraes import * # incremental conditional averaging

##################################################
### 0. Configurable parameters

## Traceset, number of traces, and S-box to attack
tracesetFilename = "traces/swaes_atmega_powertraces.npz"
sampleRange      = (950, 1150) # range of samples to attack, in the format (low, high)
N                = 2000 # number of traces to attack (less or equal to the amount of traces in the file)
offset           = 0    # trace number to start from
SboxNum          = 0    # S-box to attack, counting from 0

## Leakage model
## (these parameters correspond to function names in lracpa module)
intermediateFunction = sBoxOut                  # for CPA and LRA
leakageFunction      = leakageModelHW           # for CPA
basisFunctionsModel  = basisModelSingleBits     # for LRA

## Known key for ranking
knownKeyStr = "2B7E151628AED2A6ABF7158809CF4F3C".decode("hex") # the correct key
encrypt = True # to avoid selective commenting in the following lines below 

if encrypt: # for encryption, the first round key is as is
    knownKey = np.array(map(ord, knownKeyStr), dtype="uint8")
else:       # for decryption, need to run key expansion 
    expandedKnownKey = AES().expandKey(map(ord, knownKeyStr), 16, 16 * 11) # this returs a list
    knownKey = np.array(expandedKnownKey[176-16:177], dtype="uint8")


##################################################
### 1. Log the parameters

print "---\nAttack parameters"
print "Intermediate function   :", intermediateFunction.__name__
print "CPA leakage function    :", leakageFunction.__name__
print "LRA basis functions     :", basisFunctionsModel.__name__
print "Encryption              :", encrypt
print "S-box number            :", SboxNum
print "Known roundkey          : 0x%s" % str(bytearray(knownKey)).encode("hex")


#################################################
### 2. Load samples and data

# Readout
print "---\nLoading " + tracesetFilename
t0 = time.clock()
npzfile = np.load(tracesetFilename)
data = npzfile['data'][offset:offset + N,SboxNum] # selecting only the required byte
traces = npzfile['traces'][offset:offset + N,sampleRange[0]:sampleRange[1]]
t1 = time.clock()
timeLoad = t1 - t0

# Log traceset parameters
(numTraces, traceLength) = traces.shape
print "Number of traces loaded :", numTraces
print "Trace length            :", traceLength
print "Loading time            : %0.2f s" % timeLoad

#################################################
### 2. LRA and CPA with a fixed number of traces

print "---\nAttacks with %d traces" % numTraces
print "Running CPA...",
t0 = time.clock()
CorrTraces = cpaAES(data, traces, intermediateFunction, leakageFunction)
t1 = time.clock()
print "done in %f s" % (t1 - t0)
print "Running LRA...",
t0 = time.clock()
(R2, coefs) = lraAES(data, traces, intermediateFunction, basisFunctionsModel)
t1 = time.clock()
print "done in %f s" % (t1 - t0)
print "Normalizing LRA results...",
R2norm = normalizeR2Traces(R2)
print "done"

print "---\nAttacks with %d traces and conditional averaging" % numTraces
print "Performing conditional trace averaging...",
t0 = time.clock()
(avdata, avtraces) = conditionalAveragingAESSbox(data[0:numTraces], traces[0:numTraces])
t1 = time.clock()
print "done in %f s" % (t1 - t0)
print "Running CPA on averaged traces...",
t0 = time.clock()
CorrTracesAv = cpaAES(avdata, avtraces, intermediateFunction, leakageFunction)
t1 = time.clock()
print "done in %f s" % (t1 - t0)
print "Running LRA on averaged traces...",
t0 = time.clock()
(R2Av, coefsav) = lraAES(avdata, avtraces, intermediateFunction, basisFunctionsModel)
t1 = time.clock()
print "done in %f s" % (t1 - t0)
print "Normalizing LRA results...",
R2AvNorm = normalizeR2Traces(R2Av)
print "done"

### 3.  visualize the result, highlighting the correct trace

print "---\nPlotting..."
fig, ax = plt.subplots(3,2,sharex=True, squeeze=True)

WrongKeyRange = range(0, knownKey[SboxNum]) + range(knownKey[SboxNum] + 1, 256)

ax[0][0].plot(CorrTraces[WrongKeyRange, :].T, color = 'grey')
ax[0][0].plot(CorrTraces[knownKey[SboxNum], :], 'r')

ax[1][0].plot(R2[WrongKeyRange, :].T, color = 'grey')
ax[1][0].plot(R2[knownKey[SboxNum], :], 'r')

ax[2][0].plot(R2norm[WrongKeyRange, :].T, color = 'grey')
ax[2][0].plot(R2norm[knownKey[SboxNum], :], 'r')

ax[0][1].plot(CorrTracesAv[WrongKeyRange, :].T, color = 'grey')
ax[0][1].plot(CorrTracesAv[knownKey[SboxNum], :], 'r')

ax[1][1].plot(R2Av[WrongKeyRange, :].T, color = 'grey')
ax[1][1].plot(R2Av[knownKey[SboxNum], :], 'r')

ax[2][1].plot(R2AvNorm[WrongKeyRange, :].T, color = 'grey')
ax[2][1].plot(R2AvNorm[knownKey[SboxNum], :], 'r')

# same vertical scales for correlation and R2
ax[0][0].set_ylim(ax[0][1].get_ylim())
ax[1][0].set_ylim(ax[1][1].get_ylim())

fig.suptitle("CPA and LRA on %d traces" % numTraces)
ax[0][0].set_title('Without cond. averaging')
ax[0][1].set_title('With cond. averaging')
ax[0][0].set_ylabel('Correlation')
ax[1][0].set_ylabel('R2')
ax[2][0].set_ylabel('Normalized R2')
ax[2][0].set_xlabel('Time sample')
ax[2][1].set_xlabel('Time sample')

plt.show()