'''
The serious attack on DES round in XOR out.

The code should be self-explanatory (especially if you look into lracpa.py module).

In the plots:
- red trace is for known correct candidate
- blue trace is for the winning candidate (e.g. the one with maximum peak)
- grey traces are for all other candiadte

pysca toolbox
Version: 0.3, 2015-10-22
Started by Ilya on 2014-11-26
'''

import numpy as np
import matplotlib.pyplot as plt
import struct
import time

from desutils import * # my DES utilities
from lracpa import * # my LRA-CPA toolbox
from condaverdes import * # incremental conditional averaging


##################################################
### 0. Configurable parameters

## Traceset, number of traces, and S-box to attack
tracesetFilename = "traces/des_card8_pa.npz"
sampleRange      = (400, 450) # range of smaples to attack
N                = 10000 # number of traces to attack (less or equal to the amount of traces in the file)
offset           = 0     # trace number to start from
evolutionStep    = 500   # step for intermediate reports
SboxNum          = 1     # S-box to attack, counting from 0

## Leakage model
## (these parameters correspond to function names in lracpa module)
averagingFunction    = roundXOR_valueForAveraging # for CPA and LRA
intermediateFunction = roundXOR_targetVariable    # for CPA and LRA
leakageFunction      = leakageModelHW             # for CPA
basisFunctionsModel  = basisModelSingleBits       # for LRA

## Known key for ranking
knownKey = 0x8A7400A03230DA28
encrypt = True

# get the known key
roundKeyNum = 1
if (encrypt == False):
    roundKeyNum = 16
roundKey = computeRoundKeys(knownKey, roundKeyNum)[roundKeyNum-1]
knownKeyChunk = roundKeyChunk(roundKey, SboxNum)

##################################################
### 1. Log the parameters

print "---\nAttack parameters"
print "Averaging function      :", averagingFunction.__name__
print "Intermediate function   :", intermediateFunction.__name__
print "CPA leakage function    :", leakageFunction.__name__
print "LRA basis functions     :", basisFunctionsModel.__name__
print "Encryption              :", encrypt
print "S-box number            :", SboxNum
print "Known key               : " + format(knownKey, "#018x")
print "Known round key         : " + format(roundKey, '#014x'),
print '[',
for i in range(8):
    print format(roundKeyChunk(roundKey, i), '#04x'),
print ']'


#################################################
### 2. Load samples and data

# Readout
print "---\nLoading " + tracesetFilename
t0 = time.clock()
npzfile = np.load(tracesetFilename)
data = npzfile['data'][0:N]
traces = npzfile['traces'][0:N,sampleRange[0]:sampleRange[1]]
t1 = time.clock()
timeLoad = t1 - t0

# convert data byte arrays to integers (more convenient for DES)
print "Converting data..."
datanew = []
for i in range(0, len(data)):
    datanew.append(struct.unpack('!Q', data[i][0:8].tostring())[0])
data = datanew # old data will be garbage-collected

# Log traceset parameters
(numTraces, traceLength) = traces.shape
print "Number of traces loaded :", numTraces
print "Trace length            :", traceLength
print "Loading time            : %0.2f s" % timeLoad

#################################################
### 3. Attack with fixed amount of traces
'''
print "---\nAttack"

# perform conditional averaging
CondAver = ConditionalAveragerDes(1024, traceLength)
for i in range(N):
    CondAver.addTrace(data[i], traces[i], averagingFunction, SboxNum)
(avdata, avtraces) = CondAver.getSnapshot()

# CPA
CorrTraces = cpaDES(avdata, avtraces, intermediateFunction, SboxNum, leakageFunction)

# LRA
R2, coefs = lraDES(avdata, avtraces, intermediateFunction, SboxNum, basisFunctionsModel)

### visualize results

fig = plt.figure()

# allocate grid
axCPA = plt.subplot2grid((3, 1), (0, 0))
axLRA = plt.subplot2grid((3, 1), (1, 0))
axLRAcoefs = plt.subplot2grid((3, 1), (2, 0))

# CPA
axCPA.plot(CorrTraces.T, color = 'grey')
axCPA.plot(CorrTraces[knownKeyChunk, :], 'r')
axCPA.set_xlim([0, traceLength])

# LRA
axLRA.plot(R2.T, color = 'grey')
axLRA.plot(R2[knownKeyChunk, :], 'r')
axLRA.set_xlim([0, traceLength])

# LRA coefs
coefsKnownKey = np.array(coefs[knownKeyChunk])
axLRAcoefs.pcolormesh(coefsKnownKey[:,:-1].T)
axLRAcoefs.set_xlim([0, traceLength])

# labels
fig.suptitle("CPA and LRA on %d traces" % N)
axCPA.set_ylabel('Correlation')
axLRA.set_ylabel('R2')
axLRAcoefs.set_ylabel('Basis function (bit)')
axLRAcoefs.set_xlabel('Time sample')

plt.show()
'''
#################################################
### 4. Attack with evolving amount of traces

print "---\nAttack"

t0 = time.clock()

# initialize the incremental averager
CondAver = ConditionalAveragerDes(1024, traceLength)

# allocate arrays for storing key rank evolution
numSteps = int(np.ceil(N / np.double(evolutionStep)))
keyRankEvolutionCPA = np.zeros(numSteps)
keyRankEvolutionLRA = np.zeros(numSteps)

# the incremental loop
tracesToSkip = 20 # warm-up to avoid numerical problems for small evolution step
for i in range (0, tracesToSkip - 1):
    CondAver.addTrace(data[i], traces[i], averagingFunction, SboxNum)
for i in range(tracesToSkip - 1, N):
    CondAver.addTrace(data[i], traces[i], averagingFunction, SboxNum)

    if (((i + 1) % evolutionStep == 0) or ((i + 1) == N)):

        (avdata, avtraces) = CondAver.getSnapshot()
        
        CorrTraces = cpaDES(avdata, avtraces, intermediateFunction, SboxNum, leakageFunction)
        R2, coefs = lraDES(avdata, avtraces, intermediateFunction, SboxNum, basisFunctionsModel)
        #R2 = normalizeR2Traces(R2)

        print "---\nResults after %d traces" % (i + 1)
        print "CPA"
        CorrPeaks = np.max(np.abs(CorrTraces), axis=1) # global maximization, absolute value!
        CpaWinningCandidate = np.argmax(CorrPeaks)
        CpaWinningCandidatePeak = np.max(CorrPeaks)
        CpaCorrectCandidateRank = np.count_nonzero(CorrPeaks >= CorrPeaks[knownKeyChunk])
        CpaCorrectCandidatePeak = CorrPeaks[knownKeyChunk]
        print "Winning candidate: 0x%02x, peak magnitude %f" % (CpaWinningCandidate, CpaWinningCandidatePeak)
        print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKeyChunk, CpaCorrectCandidatePeak, CpaCorrectCandidateRank)

        print "LRA"
        R2Peaks = np.max(R2, axis=1) # global maximization
        LraWinningCandidate = np.argmax(R2Peaks)
        LraWinningCandidatePeak = np.max(R2Peaks)
        LraCorrectCandidateRank = np.count_nonzero(R2Peaks >= R2Peaks[knownKeyChunk])
        LraCorrectCandidatePeak = R2Peaks[knownKeyChunk]
        print "Winning candidate: 0x%02x, peak magnitude %f" % (LraWinningCandidate, LraWinningCandidatePeak)
        print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKeyChunk, LraCorrectCandidatePeak, LraCorrectCandidateRank)

        stepCount = int(np.floor(i / np.double(evolutionStep)))
        keyRankEvolutionCPA[stepCount] = CpaCorrectCandidateRank
        keyRankEvolutionLRA[stepCount] = LraCorrectCandidateRank

t1 = time.clock()
timeAll = t1 - t0

print "---\nCumulative timing"
print "%0.2f s" % timeAll

# save the rank evolution for later processing
#np.savez("results/keyRankEvolutionSbox%02d" % SboxNum, kreCPA=keyRankEvolutionCPA, kreLRA=keyRankEvolutionLRA, step=evolutionStep)

#################################################
### 5. Visualize results

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
if CpaWinningCandidate != knownKeyChunk:
    axCPA.plot(CorrTraces[CpaWinningCandidate, :], 'blue')
axCPA.plot(CorrTraces[knownKeyChunk, :], 'r')
axRankEvolution.plot(traceNumbers, keyRankEvolutionCPA, color = 'green')
axCPA.set_xlim([0, traceLength])

# LRA
axLRA.plot(R2.T, color = 'grey')
if LraWinningCandidate != knownKeyChunk:
    axLRA.plot(R2[LraWinningCandidate, :], 'blue')
axLRA.plot(R2[knownKeyChunk, :], 'r')
axRankEvolution.plot(traceNumbers, keyRankEvolutionLRA, color = 'magenta')
axLRA.set_xlim([0, traceLength])

# LRA coefs
coefsKnownKey = np.array(coefs[knownKeyChunk])
axLRAcoefs.pcolormesh(coefsKnownKey[:,:-1].T, cmap="jet")
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
axRankEvolution.set_xlim([traceNumbers[int(np.ceil(tracesToSkip / np.double(evolutionStep))) - 1], N])
axRankEvolution.set_ylim([0, 64])
axRankEvolution.grid(b=True, which='both', color='0.65',linestyle='-')
#axRankEvolution.ticklabel_format(style='sci', axis='x', scilimits=(0,0), useOffset=True)

# Legend for rank evolution plot
axRankEvolution.legend(['CPA', 'LRA'], loc='upper right')

plt.show()