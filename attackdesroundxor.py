'''
The serious attack on DES round in XOR out.

The code should be self-explanatory (especially if you look into lracpa.py module).

In the plots:
- red trace is for known correct candidate
- blue trace is for the winning candidate (e.g. the one with maximum peak)
- grey traces are for all other candiadte

Version: 0.2, 2015-10-20
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
N                = 4000  # number of traces to attack (less or equal to the amount of traces in the file)
offset           = 0   # trace number to start from
#evolutionStep    = 100   # step for intermediate reports
SboxNum          = 0     # S-box to attack, counting from 0

## Leakage model
## (these parameters correspond to function names in lracpa module)
averagingFunction    = roundXOR_valueForAveraging # for CPA and LRA
intermediateFunction = roundXOR_targetVariable    # for CPA and LRA
leakageFunction      = leakageModelHW             # for CPA
basisFunctionsModel  = basisModelSingleBits       # for LRA

## Known key for ranking
knownKey = 0x8A7400A03230DA28
encrypt = True

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
### 3. CPA with fixed amount of traces

print "---\nAttack" 

# get the known key
roundKey = computeRoundKeys(knownKey, 1)[0]
knownKeyChunk = roundKeyChunk(roundKey, SboxNum)

# perform conditional averaging
CondAver = ConditionalAveragerDes(1024, traceLength)
for i in range(N):
    CondAver.addTrace(data[i], traces[i], averagingFunction, SboxNum)
(avdata, avtraces) = CondAver.getSnapshot()

# CPA
CorrTraces = cpaDESwithAveraging(avdata, avtraces, intermediateFunction, SboxNum, leakageFunction)

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

