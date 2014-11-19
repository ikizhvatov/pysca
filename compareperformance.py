'''
Compare perfromance of LRA and CPA with and without conditional averaging
'''

from lracpa import *

import numpy as np
import matplotlib.pyplot as plt
import time

from aes import AES # SlowAES as AES toolbox

### 0. Prerequisites

# the correct key byte (for later metrics)
correctKey = np.array([0x2B,0x7E,0x15,0x16,0x28,0xAE,0xD2,0xA6,0xAB,0xF7,0x15,0x88,0x09,0xCF,0x4F,0x3C])

### 1. Load samples and data

# parameters
SboxNum = 2

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