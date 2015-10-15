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

#from des_block import * # the DPA contest code as DES toolbox
from desutils import * # my DES utilities
from lracpa import * # my LRA-CPA toolbox
#from condaverdes import * # incremental conditional averaging


##################################################
### 0. Configurable parameters

## Traceset, number of traces, and S-box to attack
tracesetFilename = "traces/des_card8_pa.npz"
sampleRange      = (400, 450) # range of smaples to attack
N                = 4000  # number of traces to attack (less or equal to the amount of traces in the file)
evolutionStep    = 100   # step for intermediate reports
SboxNum          = 1     # S-box to attack, counting from 0

## Leakage model
## (these parameters correspond to function names in lracpa module)
#averagingFunction    = roundXOR_valueForAveraging # for CPA and LRA
#intermediateFunction = roundXOR_targetVariable    # for CPA and LRA
leakageFunction      = leakageModelHW             # for CPA

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
traces = npzfile['traces'][0:N,sampleRange[0]:sampleRange[1]]
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
### 3. CPA with fixed amount of traces

print "---\nAttack" 

#CondAver = ConditionalAveragerDes(1024, traceLength)
#for i in range(N):
#    CondAver.addTrace(data[i], traces[i], averagingFunction, SboxNum)

#(avdata, avtraces) = CondAver.getSnapshot()
#plt.plot(avtraces.T)

CorrTraces = cpaDESwithAveraging(data, traces, roundXOR_allInOne, SboxNum, leakageFunction)

plt.plot(CorrTraces.T, color = 'grey')
plt.show()
