'''
This file is part of pysca toolbox, license is GPLv3, see https://www.gnu.org/licenses/gpl-3.0.en.html
Author: Ilya Kizhvatov
Version: 1.0, 2017-05-14

Conditional averaging for DES
'''

import numpy as np

class ConditionalAveragerDes:

    def __init__(self, numValues, traceLength):
        '''Allocate the matrix of averaged traces'''
        self.avtraces = np.zeros((numValues, traceLength))
        self.counters = np.zeros(numValues)
        print 'ConditionalAverager: initialized for %d values and trace length %d' % (numValues, traceLength)

    def addTrace(self, data, trace, dataFunction, sBoxNumber):
        '''Add a single trace with corresponding single chunk of data computed based on the given function'''

        x = dataFunction(data, sBoxNumber)

        if (self.counters[x] == 0):
            self.avtraces[x] = trace
        else:
            self.avtraces[x] = self.avtraces[x] + (trace - self.avtraces[x]) / self.counters[x]
        self.counters[x] += 1

    def getSnapshot(self):
        ''' return a snapshot of the average matrix'''
        avdataSnap = np.flatnonzero(self.counters)   # get an vector of only _observed_ values
        avtracesSnap = self.avtraces[avdataSnap]     # remove lines corresponding to non-observed values
        return avdataSnap, avtracesSnap
