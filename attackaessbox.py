'''
The serious attack on the AES S-box.

The code should be self-explanatory (especially if you look into lracpa.py module).

In the plots:
- red trace is for known correct candidate
- blue trace is for the winning candidate (e.g. the one with maximum peak)
- grey traces are for all other candiadte

Version: 0.1, 2014-11-19
Started by Ilya on 2014-11-18
'''

import numpy as np
import matplotlib.pyplot as plt
import time

from aes import AES    # interweb's SlowAES toolbox
from lracpa import *   # my LRA-CPA toolbox
from condaver import * # incremental conditional averaging


##################################################
### 0. Configurable parameters

## Traceset, number of traces, and S-box to attack
tracesetFilename = "traces/swaes_atmega_powertraces.npz"
N = 50              # number of traces to attack (less or equal to the amount o f traces in the file)
evolutionStep = 2   # step for intermediate reports
SboxNum = 1         # S-box to attack, counting from 0

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
### 3. LRA and CPA with evolving amount of traces

print "---\nAttack" 

# initialize the incremental averager
CondAver = ConditionalAverager(256, traceLength)

# allocate arrays for storing key rank evolution
numSteps = np.ceil(N / np.double(evolutionStep))
keyRankEvolutionCPA = np.zeros(numSteps)
keyRankEvolutionLRA = np.zeros(numSteps)

# the incremental loop
tracesToSkip = 20 # warm-up to avoid numerical problems
for i in range (0, tracesToSkip):
    CondAver.addTrace(data[i], traces[i])
for i in range(tracesToSkip, N):
    CondAver.addTrace(data[i], traces[i])

    if ((i % evolutionStep == 0) or (i == N-1)):

        (avdata, avtraces) = CondAver.getSnapshot()

        CorrTraces = cpaAES(avdata, avtraces, intermediateFunction, leakageFunction)
        R2 = lraAES(avdata, avtraces, intermediateFunction, basisFunctionsModel)

        print "---\nResults after %d traces" % i
        print "CPA"
        CorrPeaks = np.max(CorrTraces, axis=1) # global maximization
        CpaWinningCandidate = np.argmax(CorrPeaks)
        CpaWinningCandidatePeak = np.max(CorrPeaks)
        CpaCorrectCandidateRank = np.count_nonzero(CorrPeaks >= CorrPeaks[knownKey[SboxNum]])
        CpaCorrectCandidatePeak = CorrPeaks[knownKey[SboxNum]]
        print "Winning candidate: 0x%02x, peak magnitude %f" % (CpaWinningCandidate, CpaWinningCandidatePeak)
        print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKey[SboxNum], CpaCorrectCandidatePeak, CpaCorrectCandidateRank)

        print "LRA"
        R2Peaks = np.max(R2, axis=1) # global maximization
        LraWinningCandidate = np.argmax(R2Peaks)
        LraWinningCandidatePeak = np.max(R2Peaks)
        LraCorrectCandidateRank = np.count_nonzero(R2Peaks >= R2Peaks[knownKey[SboxNum]])
        LraCorrectCandidatePeak = R2Peaks[knownKey[SboxNum]]
        print "Winning candidate: 0x%02x, peak magnitude %f" % (LraWinningCandidate, LraWinningCandidatePeak)
        print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKey[SboxNum], LraCorrectCandidatePeak, LraCorrectCandidateRank)

        stepCount = np.ceil(i / evolutionStep)
        keyRankEvolutionCPA[stepCount] = CpaCorrectCandidateRank
        keyRankEvolutionLRA[stepCount] = LraCorrectCandidateRank


#################################################
### 4. Visualize results

print "---\nPlotting..."
fig, ax = plt.subplots(2, 2, sharex=False, squeeze=True)

# CPA
ax[0][0].plot(CorrTraces.T, color = 'grey')
ax[0][0].plot(CorrTraces[CpaWinningCandidate, :], 'blue')
ax[0][0].plot(CorrTraces[knownKey[SboxNum], :], 'r')

ax[0][1].plot(keyRankEvolutionCPA, color = 'black')

# LRA
ax[1][0].plot(R2.T, color = 'grey')
ax[1][0].plot(R2[LraWinningCandidate, :], 'blue')
ax[1][0].plot(R2[knownKey[SboxNum], :], 'r')

ax[1][1].plot(keyRankEvolutionLRA, color = 'black')

# labels
fig.suptitle("CPA and LRA on %d traces" % N)
ax[0][0].set_ylabel('Correlation')
ax[1][0].set_ylabel('R2')
ax[1][0].set_xlabel('Time sample')
ax[0][1].set_ylabel('Rank')
ax[1][1].set_ylabel('Rank')
ax[1][1].set_xlabel('Number of traces')

# Limits
ax[0][1].set_xlim([np.ceil(tracesToSkip / evolutionStep), numSteps])
ax[0][1].set_xlim([np.ceil(tracesToSkip / evolutionStep), numSteps])
ax[1][1].set_ylim([0, 256])
ax[1][1].set_ylim([0, 256])

plt.show()
