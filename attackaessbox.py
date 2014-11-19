'''
The serious attack on the AES S-box.

The code should be self-explanatory (especially if you look into lracpa.py module).

In the plots:
- red trace is for known correct candidate
- blue trace is for the winning candidate (e.g. the one with maximum peak)
- grey traces are for all other candiadte

Started by Ilya on 2014-11-18
'''

import numpy as np
import matplotlib.pyplot as plt
import time

from aes import AES   # interweb's SlowAES toolbox
from lracpa import *  # my LRA-CPA toolbox


##################################################
### 0. Configurable parameters

## Traceset, number of traces, and S-box to attack
tracesetFilename = "traces/swaes_atmega_powertraces.npz"
N = 20     # number of traces to attack (less or equal to the amount o f traces in the file)
SboxNum = 1  # S-box to attack, counting from 0

## Leakage model
## (these parameters correspond to function names in lracpa module)
intermediateFunction = sBoxOut         # for CPA and LRA
leakageFunction      = leakageModelHW  # for CPA
basisFunctionsModel  = basisModel9     # for LRA

## Known key for ranking
knownKeyStr = "2B7E151628AED2A6ABF7158809CF4F3C".decode("hex") # the correct key
encrypt = True # to avoid selective commenting in the followign lines below 

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
data = npzfile['data'][0:N,SboxNum] # selecting only the required byte
traces = npzfile['traces'][0:N,:]   # can select a sub-range of samples here TODO: parametrize!
t1 = time.clock()
timeLoad = t1 - t0

# Log traceset parameters
(numTraces, traceLength) = traces.shape
print "Number of traces loaded :", numTraces
print "Trace length            :", traceLength
print "Loading time            : %0.2f s" % timeLoad


#################################################
### 3. LRA and CPA with a fixed number of traces

print "---\nAttack" 

print "Performing conditional trace averaging...",
t0 = time.clock()
(avdata, avtraces) = conditionalAveragingAESSbox(data, traces)
#plt.plot(avtraces.T)
#plt.show()
t1 = time.clock()
print "done in %0.2f s" % (t1 - t0)

print "Running CPA on averaged traces...",
t0 = time.clock()
CorrTraces = cpaAES(avdata, avtraces, intermediateFunction, leakageFunction)
t1 = time.clock()
print "done in %0.2f s" % (t1 - t0)

print "Running LRA on averaged traces...",
t0 = time.clock()
R2 = lraAES(avdata, avtraces, intermediateFunction, basisFunctionsModel)
t1 = time.clock()
print "done in %0.2f s" % (t1 - t0)


#################################################
### 4. Analyze results

print "---\nResults CPA"
CorrPeaks = np.max(CorrTraces, axis=1) # global maximization
CpaWinningCandidate = np.argmax(CorrPeaks)
CpaWinningCandidatePeak = np.max(CorrPeaks)
CpaCorrectCandidateRank = np.count_nonzero(CorrPeaks >= CorrPeaks[knownKey[SboxNum]])
CpaCorrectCandidatePeak = CorrPeaks[knownKey[SboxNum]]
print "Winning candidate: 0x%02x, peak magnitude %f" % (CpaWinningCandidate, CpaWinningCandidatePeak)
print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKey[SboxNum], CpaCorrectCandidatePeak, CpaCorrectCandidateRank)

print "---\nResults LRA"
R2Peaks = np.max(R2, axis=1) # global maximization
LraWinningCandidate = np.argmax(R2Peaks)
LraWinningCandidatePeak = np.max(R2Peaks)
LraCorrectCandidateRank = np.count_nonzero(R2Peaks >= R2Peaks[knownKey[SboxNum]])
LraCorrectCandidatePeak = R2Peaks[knownKey[SboxNum]]
print "Winning candidate: 0x%02x, peak magnitude %f" % (LraWinningCandidate, LraWinningCandidatePeak)
print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKey[SboxNum], LraCorrectCandidatePeak, LraCorrectCandidateRank)


#################################################
### 5. Visualize results

print "---\nPlotting..."
fig, ax = plt.subplots(2, 1, sharex=True, squeeze=True)

# CPA
ax[0].plot(CorrTraces.T, color = 'grey')
ax[0].plot(CorrTraces[CpaWinningCandidate, :], 'blue')
ax[0].plot(CorrTraces[knownKey[SboxNum], :], 'r')

# LRA
ax[1].plot(R2.T, color = 'grey')
ax[1].plot(R2[LraWinningCandidate, :], 'blue')
ax[1].plot(R2[knownKey[SboxNum], :], 'r')

# labels
fig.suptitle("CPA and LRA on %d traces" % N)
ax[0].set_ylabel('Correlation')
ax[1].set_ylabel('R2')
ax[1].set_xlabel('Time sample')

plt.show()
