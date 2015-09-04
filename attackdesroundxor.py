'''
The serious attack on DES round in XOR out.

The code should be self-explanatory (especially if you look into lracpa.py module).

In the plots:
- red trace is for known correct candidate
- blue trace is for the winning candidate (e.g. the one with maximum peak)
- grey traces are for all other candiadte

Started by Ilya on 2014-11-26
'''

import numpy as np
import matplotlib.pyplot as plt
import struct
import time

from des_block import *  # the DPA contest code as DES toolbox
from desutils  import *  # my DES utilities
from lracpa    import *  # my LRA-CPA toolbox
from condaver  import *  # incremental conditional averaging


##################################################
### 0. Configurable parameters

## Traceset, number of traces, and S-box to attack
tracesetFilename = "traces/des_card8_pa.npz"
N                = 4000  # number of traces to attack (less or equal to the amount o f traces in the file)
evolutionStep    = 100   # step for intermediate reports
SboxNum          = 0     # S-box to attack, counting from 0

## Leakage model
## (these parameters correspond to function names in lracpa module)
averagingFunction    = roundXOR_valueForAveraging # for CPA and LRA
intermediateFunction = roundXOR_targetVariable    # for CPA and LRA
leakageFunction      = leakageModelHW             # for CPA
basisFunctionsModel  = basisModelSingleBits       # for LRA

## Known key for ranking
knownKeyStr = "8A7400A03230DA28".decode("hex") # the correct key
encrypt = True # to avoid selective commenting in the following lines below 

# TODO for DES: compute the 1st round key
#if encrypt: # for encryption, the first round key is as is
#    knownKey = np.array(map(ord, knownKeyStr), dtype="uint8")
#else:       # for decryption, need to run key expansion 
#    expandedKnownKey = AES().expandKey(map(ord, knownKeyStr), 16, 16 * 11) # this returs a list
#    knownKey = np.array(expandedKnownKey[176-16:177], dtype="uint8")


##################################################
### 1. Log the parameters

print "---\nAttack parameters"
print "Averaging function      :", averagingFunction.__name__
print "Intermediate function   :", intermediateFunction.__name__
print "CPA leakage function    :", leakageFunction.__name__
print "LRA basis functions     :", basisFunctionsModel.__name__
print "Encryption              :", encrypt
print "S-box number            :", SboxNum
#print "Known roundkey          : 0x%s" % str(bytearray(knownKey)).encode("hex") TODO


#################################################
### 2. Load samples and data

# Readout
print "---\nLoading " + tracesetFilename
t0 = time.clock()
npzfile = np.load(tracesetFilename)
data = npzfile['data'][0:N] # selecting only the required byte
traces = npzfile['traces'][0:N,:]   # can select a sub-range of samples here TODO: parametrize!
t1 = time.clock()
timeLoad = t1 - t0

# convert data byte arrays to integers (more convenient for DES)
print "Gathering bytes to Python long integers (no numpy uint64 as numpy shifts do not support it!)"
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
### 3. LRA and CPA with evolving amount of traces

print "---\nAttack" 

# initialize the incremental averager
# TODO: 1024 depends on the averagingFunction, needs to be parametrized
CondAver = ConditionalAveragerDes(1024, traceLength)

# allocate arrays for storing key rank evolution
numSteps = np.ceil(N / np.double(evolutionStep))
keyRankEvolutionCPA = np.zeros(numSteps)
keyRankEvolutionLRA = np.zeros(numSteps)

# the incremental loop
tracesToSkip = 50 # warm-up to avoid numerical problems
for i in range (0, tracesToSkip - 1):
    CondAver.addTrace(data[i], traces[i], averagingFunction, SboxNum)
for i in range(tracesToSkip - 1, N):
    CondAver.addTrace(data[i], traces[i], averagingFunction, SboxNum)

    if ((i % evolutionStep == 0) or (i == N-1)):

        (avdata, avtraces) = CondAver.getSnapshot()

        CorrTraces = cpaDESwithAveraging(avdata, avtraces, intermediateFunction, SboxNum, leakageFunction)
        R2 = np.zeros((traceLength,64))#lraDES(avdata, avtraces, intermediateFunction, SboxNum, basisFunctionsModel)

        print "---\nResults after %d traces" % i
        print "CPA"
        CorrPeaks = np.max(CorrTraces, axis=1) # global maximization
        CpaWinningCandidate = np.argmax(CorrPeaks)
        CpaWinningCandidatePeak = np.max(CorrPeaks)
        #CpaCorrectCandidateRank = np.count_nonzero(CorrPeaks >= CorrPeaks[knownKey[SboxNum]])
        #CpaCorrectCandidatePeak = CorrPeaks[knownKey[SboxNum]]
        print "Winning candidate: 0x%02x, peak magnitude %f" % (CpaWinningCandidate, CpaWinningCandidatePeak)
        #print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKey[SboxNum], CpaCorrectCandidatePeak, CpaCorrectCandidateRank)

        print "LRA"
        R2Peaks = np.max(R2, axis=1) # global maximization
        LraWinningCandidate = np.argmax(R2Peaks)
        LraWinningCandidatePeak = np.max(R2Peaks)
        #LraCorrectCandidateRank = np.count_nonzero(R2Peaks >= R2Peaks[knownKey[SboxNum]])
        #LraCorrectCandidatePeak = R2Peaks[knownKey[SboxNum]]
        print "Winning candidate: 0x%02x, peak magnitude %f" % (LraWinningCandidate, LraWinningCandidatePeak)
        #print "Correct candidate: 0x%02x, peak magnitude %f, rank %d" % (knownKey[SboxNum], LraCorrectCandidatePeak, LraCorrectCandidateRank)

        stepCount = np.ceil(i / evolutionStep)
        #keyRankEvolutionCPA[stepCount] = CpaCorrectCandidateRank
        #keyRankEvolutionLRA[stepCount] = LraCorrectCandidateRank


#################################################
### 4. Visualize results

print "---\nPlotting..."

fig = plt.figure()

# allocate grid
axCPA = plt.subplot2grid((2, 2), (0, 0))
axLRA = plt.subplot2grid((2, 2), (1, 0))
axRankEvolution = plt.subplot2grid((2, 2), (0, 1), rowspan = 2)

# CPA
axCPA.plot(CorrTraces.T, color = 'grey')
#if CpaWinningCandidate == knownKey[SboxNum]:
axCPA.plot(CorrTraces[CpaWinningCandidate, :], 'blue')
#axCPA.plot(CorrTraces[knownKey[SboxNum], :], 'r')
axRankEvolution.plot(keyRankEvolutionCPA, color = 'green')

# LRA
#axLRA.plot(R2.T, color = 'grey')
#if LraWinningCandidate == knownKey[SboxNum]:
#    axLRA.plot(R2[LraWinningCandidate, :], 'blue')
#axLRA.plot(R2[knownKey[SboxNum], :], 'r')
#axRankEvolution.plot(keyRankEvolutionLRA, color = 'magenta')

# labels
fig.suptitle("CPA and LRA on %d traces" % N)
axCPA.set_ylabel('Correlation')
axLRA.set_ylabel('R2')
axLRA.set_xlabel('Time sample')
axRankEvolution.set_ylabel('Correct key candidate rank')
axRankEvolution.set_xlabel('Number of traces')
axRankEvolution.set_title('Correct key rank evolution (global maximisation)')

# Limits and tick labels for key rand evolution plot
# TODO fix label value offset (should be -1 to what it is now)
axRankEvolution.set_xlim([np.ceil(tracesToSkip / np.double(evolutionStep)), numSteps - 1])
axRankEvolution.set_ylim([0, 64])
xTickLocations = axRankEvolution.get_xticks()
xTickLabels = map(lambda x: "%g" % x, xTickLocations * evolutionStep)
axRankEvolution.set_xticklabels(xTickLabels)

# Legend for rank evolution plot
axRankEvolution.legend(['CPA', 'LRA'], loc='upper right')

plt.show()
