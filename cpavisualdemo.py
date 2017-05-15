'''
This file is part of pysca toolbox, license is GPLv3, see https://www.gnu.org/licenses/gpl-3.0.en.html
Author: Ilya Kizhvatov
Version: 1.0, 2017-05-14

CPA visual demo with live evolution, with AES-128 traces from ATmega16
'''

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # enforcing the backend, otherwise fails on Mac OS X with its default "macosx" backend
import matplotlib.pyplot as plt
import matplotlib.patches as pch

from lracpa import correlationTraces # my LRA-CPA toolbox

##############################################################################
# Parameters

numByte = 0 # byte to attack, from 0 to 15

# range for the nubmer of traces to attack (there are 2000 traces)
startTraces = 5 # less than 5 leads to numerical problems in correlation computation
stopTraces = 100
stepTraces = 1

#  range of samples to attack (traces are 2800 samples long)
sampleRange = range(950,1100)


##############################################################################
# Prerequisites

# the correct key (for later metrics, if any)
correctKey = np.array([0x2B,0x7E,0x15,0x16,0x28,0xAE,0xD2,0xA6,0xAB,0xF7,0x15,0x88,0x09,0xCF,0x4F,0x3C])

# hamming weight leakage model (pure)
byteHammingWeight = np.load('data/bytehammingweight.npy') # HW table
def leakageModel_HW(x):
    return byteHammingWeight[x]


##############################################################################
# Preload and precompute

# load AES S-Box
sbox = np.load('data/aessbox.npy')

# load samples and data
npzfile = np.load('traces/swaes_atmega_power.npz')
data = npzfile['data']
traces = npzfile['traces'][:,sampleRange]

# output traceset parameters
numTraces = traces.shape[0]
traceLength = traces.shape[1]
print "Number of traces: ", numTraces
print "Trace length: ", traceLength

# compute intermediate variable hypotheses for all the key candidates
k = np.arange(0,256, dtype='uint8')
H = np.zeros((256, len(data)), dtype='uint8')
for i in range(256):
    H[i,:] = sbox[data[:, numByte] ^ k[i]]

# compute leakage hypotheses for every 
HL = np.array(map(leakageModel_HW, H)).T # leakage model here can be changed


##############################################################################
# CPA attack (interleaved with incremental plotting, so a bit of a mess)

### Graphics stuff
# allocate a line object for every correlation trace and evolution trace
hc = []
he = []
fig, (ax1,ax2) = plt.subplots(1,2, sharey=True)
for i in range(256):
    hc.append(ax2.plot([],[], linewidth=2, color='grey')[0])
    he.append(ax1.plot([],[], linewidth=2, color='grey')[0])
# put the text label for showing the winning key candidate    
ht = plt.text(10, -0.95, "", fontsize=18, fontweight='bold')
ax1.set_ylabel('Correlation', fontsize=18)
ax1.set_xlabel('Number of traces', fontsize=18)
ax2.set_xlabel('Time sample', fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=18)
ax2.tick_params(axis='both', which='major', labelsize=18)
# show the figure and save the background
ax2.axis([0, traceLength, -1, 1])
ax1.axis([0, stopTraces, -1, 1])
fig.show()
fig.canvas.draw()
for i in range(256):
    hc[i].set_xdata(range(1,traceLength+1))

hp = pch.ConnectionPatch(xyA=(0,0), xyB =(0,0), coordsA='data', axesA=ax2, axesB=ax1, color='black', linestyle='dashed')
ax2.add_artist(hp)

# main loop; attack and graphics are inevitably interleaved
corrTraces = np.empty([256, traceLength]);
guessedKeyBytePrev = 0
for n in range(startTraces, stopTraces, stepTraces):

    # compute correlation traces, update them in the plot
    corrTraces = correlationTraces(traces[0:n], HL[0:n])
    for i in range(0, 256):
        hc[i].set_ydata(corrTraces[i])
        ax1.draw_artist(hc[i])

    # determine the peaks and the most probable key candidate
    corrPeaks = np.max(corrTraces, axis=1)
    guessedKeyByte = np.argmax(corrPeaks)

    # update highlighting to the current winning key candidate
    hc[guessedKeyBytePrev].set_color('grey')
    hc[guessedKeyByte].set_color('red')
    he[guessedKeyBytePrev].set_color('grey')
    he[guessedKeyByte].set_color('red')
    guessedKeyBytePrev = guessedKeyByte

    for i in range(0, 256):
        he[i].set_xdata(np.append(he[i].get_xdata(), n))
        he[i].set_ydata(np.append(he[i].get_ydata(), corrPeaks[i]))
        ax1.draw_artist(he[i])
    
    # update the label
    label = "Key 0x%02x, peak %0.4f" % (guessedKeyByte, np.max(corrPeaks))
    ht.set_text(label)
    ax2.draw_artist(ht)

    # show connecting line between the peak and the evolution plot
    hp.xy2 = (n, corrPeaks[guessedKeyByte])
    hp.xy1 = (np.argmax(corrTraces[guessedKeyByte]), corrPeaks[guessedKeyByte])
    ax2.draw_artist(hp)
    
    # refresh the plot
    # NOTE: no need to use blitting here or otherwise optimize
    #       we even do not want to be very fast!
    fig.canvas.draw()
    fig.canvas.flush_events() # crucial: prevents freezing!

plt.show() # to keep the figure open
