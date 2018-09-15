'''
This file is part of pysca toolbox, license is GPLv3, see https://www.gnu.org/licenses/gpl-3.0.en.html
Author: Ilya Kizhvatov
Version: 1.0, 2017-05-14

Convert Inspector traceset into numpy array and save to npz.
Reads the entire traceset into memory, so cannot deal with huge tracesets.

External packages required:
- numpy
- Trace.py

Versions of this file:
v0.2 2014-10-31 Ilya: data conversion to unit64 made optional; some refactoring
v0.1 2013-11-12 Ilya: initial
'''

import argparse
import numpy as np
import struct
import Trace as trs

# Determine coding of samples in .trs traceset and return it in numpy format
def determineTrsSampleCoding(ts):
    if ts._sampleCoding == ts.CodingByte:
        samplesDataType = "int8"
    elif ts._sampleCoding == ts.CodingShort:
        samplesDataType = "int16"
    elif ts._sampleCoding == ts.CodingInt:
        samplesDataType = "int32"
    elif ts._sampleCoding == ts.CodingFloat:
        samplesDataType = "float32"
    else:
        samplesDataType = None
    return samplesDataType
    
# Print main metadata of the .trs traceset
def printTrsMetadata(ts, samplesDataType):
    print("Number of traces:\t%d" % ts._numberOfTraces)
    print("Samples per trace:\t%d" % ts._numberOfSamplesPerTrace)
    print("Samples datatype:\t%s" % samplesDataType)
    print("Data bytes:\t\t%d" % ts._dataSpace)
    print("Trace block size:\t%d bytes" % ts._traceBlockSpace)
    print("Header size:\t\t%d bytes" % ts._traceBlockOffset)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert Inspector 4 traceset into numpy array')
    parser.add_argument('-c', '--convertdata',  action='store_true', help='convert data from byte array to uint64 chunks')
    parser.add_argument('filename', help='traceset file name without trs extension')
    args = parser.parse_args()

    ts = trs.TraceSet()
    ts.open(args.filename + ".trs")
    samplesDataType = determineTrsSampleCoding(ts)
    printTrsMetadata(ts, samplesDataType)

    # read out the traces
    print("Preallocating arrays")
    traces = np.empty(shape=(ts._numberOfTraces, ts._numberOfSamplesPerTrace), dtype = samplesDataType)
    data = np.empty(shape=(ts._numberOfTraces, ts._dataSpace), dtype = "uint8")
    print("Populating arrays")
    for i in range(ts._numberOfTraces):
        t = ts.getTrace(i)
        traces[i, :] = np.array(t._samples, dtype = samplesDataType)
        data[i, :] = np.array(t._data, dtype = "uint8")

    if args.convertdata:
        print("Gathering bytes to uint64's")
        wordBytes = 8
        structCommand = '!Q'
        numberOfWordsNecessary = (len(data[0]) + wordBytes - 1)//(wordBytes)
        datanew = np.empty((len(data),numberOfWordsNecessary), dtype='uint64')
        for i in range(0, len(data)):
            if(len(data[i]) < numberOfWordsNecessary*wordBytes):
                tempData = np.concatenate((np.zeros(numberOfWordsNecessary*wordBytes - len(data[i]), dtype="uint8"), data[i]))
            else:
                tempData = data[i]
            for j in range(numberOfWordsNecessary):
                datanew[i][j] = struct.unpack(structCommand, tempData[(j*wordBytes):((j+1)*wordBytes)].tostring())[0]
        data = datanew # old data will be garbage-collected

    print("Saving file")
    np.savez(args.filename, traces=traces, data=data)
    print("Done")