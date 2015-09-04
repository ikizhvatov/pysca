'''
The serious attack on the AES S-box.

The code should be self-explanatory (especially if you look into lracpa.py module).

In the plots:
- red trace is for known correct candidate
- blue trace is for the winning candidate (e.g. the one with maximum peak)
- grey traces are for all other candiadte

Version: 0.2, 2015-09-04
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
N                = 100 # number of traces to attack (less or equal to the amount o f traces in the file)
offset           = 0   # trace number to start from
evolutionStep    = 10  # step for intermediate reports
SboxNum          = 3   # S-box to attack, counting from 0

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
traces = npzfile['traces'][offset:offset + N,800:1500]   # can select a sub-range of samples here TODO: parametrize!
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

t0 = time.clock()

# initialize the incremental averager
CondAver = ConditionalAverager(256, traceLength)

# allocate arrays for storing key rank evolution
numSteps = np.ceil(N / np.double(evolutionStep))
keyRankEvolutionCPA = np.zeros(numSteps)
keyRankEvolutionLRA = np.zeros(numSteps)

# the incremental loop
tracesToSkip = 20 # warm-up to avoid numerical problems for small evolution step
for i in range (0, tracesToSkip - 1):
    CondAver.addTrace(data[i], traces[i])
for i in range(tracesToSkip - 1, N):
    CondAver.addTrace(data[i], traces[i])

    if (((i + 1) % evolutionStep == 0) or ((i + 1) == N)):

        (avdata, avtraces) = CondAver.getSnapshot()
        
        CorrTraces = cpaAES(avdata, avtraces, intermediateFunction, leakageFunction)
        R2, coefs = lraAES(avdata, avtraces, intermediateFunction, basisFunctionsModel)

        print "---\nResults after %d traces" % (i + 1)
        print "CPA"
        CorrPeaks = np.max(np.abs(CorrTraces), axis=1) # global maximization, absolute value!
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

        stepCount = np.floor(i / np.double(evolutionStep))
        keyRankEvolutionCPA[stepCount] = CpaCorrectCandidateRank
        keyRankEvolutionLRA[stepCount] = LraCorrectCandidateRank

t1 = time.clock()
timeAll = t1 - t0

#################################################
### 4. Visualize results

# save the rank evolution for later processing
#np.savez("results/keyRankEvolutionSbox%02d" % SboxNum, kreCPA=keyRankEvolutionCPA, kreLRA=keyRankEvolutionLRA, step=evolutionStep)

print "---\nCumulative timing"
print "%0.2f s" % timeAll

print "---\nPlotting..."

fig = plt.figure()

# allocate grid
axCPA = plt.subplot2grid((3, 2), (0, 0))
axLRA = plt.subplot2grid((3, 2), (1, 0))
axLRAcoefs = plt.subplot2grid((3, 2), (2, 0))
axRankEvolution = plt.subplot2grid((2, 2), (0, 1), rowspan = 3)

# compute trace nubmers for x axis (TODO: move into block 3)
traceNumbers = np.arange(evolutionStep, N + 1, evolutionStep)

# CPA
axCPA.plot(CorrTraces.T, color = 'grey')
if CpaWinningCandidate != knownKey[SboxNum]:
    axCPA.plot(CorrTraces[CpaWinningCandidate, :], 'blue')
axCPA.plot(CorrTraces[knownKey[SboxNum], :], 'r')
axRankEvolution.plot(traceNumbers, keyRankEvolutionCPA, color = 'green')
axCPA.set_xlim([0, traceLength])

# LRA
axLRA.plot(R2.T, color = 'grey')
if LraWinningCandidate != knownKey[SboxNum]:
    axLRA.plot(R2[LraWinningCandidate, :], 'blue')
axLRA.plot(R2[knownKey[SboxNum], :], 'r')
axRankEvolution.plot(traceNumbers, keyRankEvolutionLRA, color = 'magenta')
axLRA.set_xlim([0, traceLength])

# LRA coefs
coefsKnownKey = np.array(coefs[knownKey[SboxNum]])
axLRAcoefs.pcolormesh(coefsKnownKey[:,:-1].T)
axLRAcoefs.set_xlim([0, traceLength])

# labels
fig.suptitle("CPA and LRA on %d traces" % N)
axCPA.set_ylabel('Correlation')
axLRA.set_ylabel('R2')
axLRAcoefs.set_ylabel('Basis function (bit)')
axLRAcoefs.set_xlabel('Time sample')
axRankEvolution.set_ylabel('Correct key candidate rank')
axRankEvolution.set_xlabel('Number of traces')
axRankEvolution.set_title('Correct key rank evolution (global maximisation)')

# Limits and tick labels for key rand evolution plot
axRankEvolution.set_xlim([traceNumbers[np.ceil(tracesToSkip / np.double(evolutionStep)) - 1], N])
axRankEvolution.set_ylim([0, 256])
axRankEvolution.grid(b=True, which='both', color='0.65',linestyle='-')
#axRankEvolution.ticklabel_format(style='sci', axis='x', scilimits=(0,0), useOffset=True)

# Legend for rank evolution plot
axRankEvolution.legend(['CPA', 'LRA'], loc='upper right')

plt.show()
