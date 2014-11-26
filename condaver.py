'''
Conditional averaging

Very rough so far, and only for AES

TODO: 
- iterator (therefore loop inside)
- automatic readout from a file

Started by Ilya on 2014-11-20
'''

import numpy as np

class ConditionalAverager:

    def __init__(self, numValues, traceLength):
        '''Allocate the matrix of averaged traces'''
        self.avtraces = np.zeros((numValues, traceLength))
        self.counters = np.zeros(numValues)
        print 'ConditionalAverager: initialized for %d values and trace length %d' % (numValues, traceLength)

    def addTrace(self, data, trace):
        '''Add a single trace with corresponding single chunk of data'''
        if (self.counters[data] == 0):
            self.avtraces[data] = trace
        else:
            self.avtraces[data] = self.avtraces[data] + (trace - self.avtraces[data]) / self.counters[data]
        self.counters[data] += 1

    def getSnapshot(self):
        ''' return a snapshot of the average matrix'''
        avdataSnap = np.flatnonzero(self.counters)   # get an vector of only _observed_ values
        avtracesSnap = self.avtraces[avdataSnap] # remove lines corresponding to non-observed values
        return avdataSnap, avtracesSnap

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
