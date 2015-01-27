'''
DES transformations required for DPA with conditional averaging

Uses minor chunks of code from pyDES-2.0.1 and DPA contest v1 DES example.

TODO: rewrite in Cython or in C using cyclic shifts and other natural bitwise
      operations; look at DES implementation in libtomcrypt as an example.

Started by Ilya on 2014-11-25
'''

import numpy as np


##############################################################################
# Core functionality

def permuteBits(x, permutation):
    ''' Permutes bits of x given a permutation table. Implemented as in
        Inspector 4 code. Assumes that permutation table is 0-offset. '''
    result = 0L
    for i in range(0, len(permutation)):
        result = ((result << 1) |
                  ((x >> (len(permutation) - 1 - permutation[i])) & 1))
    return result

# This is a 64-bit permutaion, not a lookup table
# copied from pyDES-2.0.1
InitialPermutation = [
    57, 49, 41, 33, 25, 17,  9, 1,
    59, 51, 43, 35, 27, 19, 11, 3,
    61, 53, 45, 37, 29, 21, 13, 5,
    63, 55, 47, 39, 31, 23, 15, 7,
    56, 48, 40, 32, 24, 16,  8, 0,
    58, 50, 42, 34, 26, 18, 10, 2,
    60, 52, 44, 36, 28, 20, 12, 4,
    62, 54, 46, 38, 30, 22, 14, 6
    ]

# This is a 6-to-4 bits lookup table
# copied from pyDES-2.0.1
SBoxLUT = [
    # S1
    np.array([14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7,
     0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8,
     4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0,
     15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13], dtype='uint8'),

    # S2
    np.array([15, 1, 8, 14, 6, 11, 3, 4, 9, 7, 2, 13, 12, 0, 5, 10,
     3, 13, 4, 7, 15, 2, 8, 14, 12, 0, 1, 10, 6, 9, 11, 5,
     0, 14, 7, 11, 10, 4, 13, 1, 5, 8, 12, 6, 9, 3, 2, 15,
     13, 8, 10, 1, 3, 15, 4, 2, 11, 6, 7, 12, 0, 5, 14, 9], dtype='uint8'),

    # S3
    np.array([10, 0, 9, 14, 6, 3, 15, 5, 1, 13, 12, 7, 11, 4, 2, 8,
     13, 7, 0, 9, 3, 4, 6, 10, 2, 8, 5, 14, 12, 11, 15, 1,
     13, 6, 4, 9, 8, 15, 3, 0, 11, 1, 2, 12, 5, 10, 14, 7,
     1, 10, 13, 0, 6, 9, 8, 7, 4, 15, 14, 3, 11, 5, 2, 12], dtype='uint8'),

    # S4
    np.array([7, 13, 14, 3, 0, 6, 9, 10, 1, 2, 8, 5, 11, 12, 4, 15,
     13, 8, 11, 5, 6, 15, 0, 3, 4, 7, 2, 12, 1, 10, 14, 9,
     10, 6, 9, 0, 12, 11, 7, 13, 15, 1, 3, 14, 5, 2, 8, 4,
     3, 15, 0, 6, 10, 1, 13, 8, 9, 4, 5, 11, 12, 7, 2, 14], dtype='uint8'),

    # S5
    np.array([2, 12, 4, 1, 7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9,
     14, 11, 2, 12, 4, 7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6,
     4, 2, 1, 11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14,
     11, 8, 12, 7, 1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3], dtype='uint8'),

    # S6
    np.array([12, 1, 10, 15, 9, 2, 6, 8, 0, 13, 3, 4, 14, 7, 5, 11,
     10, 15, 4, 2, 7, 12, 9, 5, 6, 1, 13, 14, 0, 11, 3, 8,
     9, 14, 15, 5, 2, 8, 12, 3, 7, 0, 4, 10, 1, 13, 11, 6,
     4, 3, 2, 12, 9, 5, 15, 10, 11, 14, 1, 7, 6, 0, 8, 13], dtype='uint8'),

    # S7
    np.array([4, 11, 2, 14, 15, 0, 8, 13, 3, 12, 9, 7, 5, 10, 6, 1,
     13, 0, 11, 7, 4, 9, 1, 10, 14, 3, 5, 12, 2, 15, 8, 6,
     1, 4, 11, 13, 12, 3, 7, 14, 10, 15, 6, 8, 0, 5, 9, 2,
     6, 11, 13, 8, 1, 4, 10, 7, 9, 5, 0, 15, 14, 2, 3, 12], dtype='uint8'),

    # S8
    np.array([13, 2, 8, 4, 6, 15, 11, 1, 10, 9, 3, 14, 5, 0, 12, 7,
     1, 15, 13, 8, 10, 3, 7, 4, 12, 5, 6, 11, 0, 14, 9, 2,
     7, 11, 4, 1, 9, 12, 14, 2, 0, 6, 10, 13, 15, 3, 5, 8,
     2, 1, 14, 7, 4, 10, 8, 13, 15, 12, 9, 0, 3, 5, 6, 11], dtype='uint8')
]


''' Extension permutation, as a list of function per S-box.
    x is supposed to be a 32-bit wide integer. If not, function will work
    incorrectly.
    Calling example: ExtensionPermutationPerSbox[4](x) - retrive the part
    of E(x) corresponding to S-box 4'''
ExpansionPerSbox = {
    0 : lambda x : ((x >> 27) | (x << 5)) & 0x3f,
    1 : lambda x : (x >> 23) & 0x3f,
    2 : lambda x : (x >> 19) & 0x3f,
    3 : lambda x : (x >> 15) & 0x3f,
    4 : lambda x : (x >> 11) & 0x3f,
    5 : lambda x : (x >> 7) & 0x3f,
    6 : lambda x : (x >> 3) & 0x3f,
    7 : lambda x : ((x << 1) | (x >> 31)) & 0x3f
}

'''Inverse P permutation as a list of bit-gathering functions per S-box.
   Generated using a helper function below.'''
InversePermutationPerSbox = {
    0 : lambda x : ((x >> 20) & 8) | ((x >> 13) & 4) | ((x >>  8) & 2) | ((x >>  1) & 1),
    1 : lambda x : ((x >> 16) & 8) | ((x >>  2) & 4) | ((x >> 29) & 2) | ((x >> 14) & 1),
    2 : lambda x : ((x >>  5) & 8) | ((x >> 14) & 4) | ((x >>  1) & 2) | ((x >> 26) & 1),
    3 : lambda x : ((x >>  3) & 8) | ((x >> 10) & 4) | ((x >> 21) & 2) | ((x >> 31) & 1),
    4 : lambda x : ((x >> 21) & 8) | ((x >> 16) & 4) | ((x >>  6) & 2) | ((x >> 29) & 1),
    5 : lambda x : ((x >> 25) & 8) | ((x >>  1) & 4) | ((x >> 20) & 2) | ((x >> 13) & 1),
    6 : lambda x : ((x <<  3) & 8) | ((x >> 18) & 4) | ((x >>  9) & 2) | ((x >> 25) & 1),
    7 : lambda x : ((x >> 24) & 8) | ((x >>  3) & 4) | ((x >> 16) & 2) | ((x >> 11) & 1)
}


##############################################################################
# Tandem of functions for round in xor out intermediate. First functions
# obtains value for conditional averaging. Second function obtains the value
# of the target variable from that

def roundXOR_valueForAveraging(input, sBoxNumber):
    ''' Compute the value for conditional averaging from input, for a given
        S-box number '''

    # prepare the first round input halves
    permutedInput = permuteBits(input, InitialPermutation)
    rightHalf = permutedInput & 0xFFFFFFFF
    leftHalf = permutedInput >> 32

    # 1. get 6-bit first part: the input chunks per S-box based on the structre value
    # 2. get 4-bit second part: the bits from XOR of left and right input
    # 3. concatenate to 10-bit value
    # TODO: consider a structure instead of concatenation
    a = ExpansionPerSbox[sBoxNumber](rightHalf)
    b = InversePermutationPerSbox[sBoxNumber](rightHalf ^ leftHalf)
    r = (a << 4) | b 

    return r

def roundXOR_targetVariable(averagingValue, keyChunk, sBoxNumber):
    ''' Compute the intermediate variable the value used for from key chunk,
        for a given S-box number '''

    # unpack the value
    x = (averagingValue >> 4) & 0x3f
    y = averagingValue & 0xf

    # compute the intermediate value
    SBoxIn = x ^ keyChunk
    SBoxOut = SBoxLUT[sBoxNumber][SBoxIn]
    RoundInXorOutPerSBox = SBoxOut ^ y

    return RoundInXorOutPerSBox



##############################################################################
# Helper functions and tests

def generateInversePermutationPerSbox():
    ''' Helper used to generate the shifts. In the output, negative values should be manually replaced by a left shift! '''
    initialPositionsPerSbox = [
        np.array([8, 16, 22, 30]),
        np.array([12, 27, 1, 17]),
        np.array([23, 15, 29, 5]),
        np.array([25, 19, 9, 0]),
        np.array([7, 13, 24, 2]),
        np.array([3, 28, 10, 18]),
        np.array([31, 11, 21, 6]),
        np.array([4, 26, 14, 20]),
        ]
    finalPositions = np.array([28, 29, 30, 31])

    for group in initialPositionsPerSbox:
        shifts = finalPositions - group
        print "((x >> %d) & 8) | ((x >> %d) & 4) | ((x >> %d) & 2) | ((x >> %d) & 1)" % (shifts[0], shifts[1], shifts[2], shifts[3])

def testDesUtilities():
    ''' Unit test for DES utilities '''
    # TODO check against a test vector generated with an existing
    #      DES implmentation

    Input = 0xA76DB873C63FE078
    KeyChunk = 0x2B
    print "Input:", hex(Input)
    print "Key chunk:", hex(KeyChunk)

    print "Values for averaging and target variables:"
    for i in range(0, 8):
        r = roundXOR_valueForAveraging(Input, i)
        t = roundXOR_targetVariable(r, KeyChunk, i)
        print "0x%04x, 0x%04x" % (r, t)

    print "Initial permutation"
    print hex(permuteBits(Input, InitialPermutation))

    RightHalf = Input & 0xFFFFFFFF

    print "Expansion per S-box:",
    print hex(ExpansionPerSbox[0](RightHalf)),
    print hex(ExpansionPerSbox[1](RightHalf)),
    print hex(ExpansionPerSbox[2](RightHalf)),
    print hex(ExpansionPerSbox[3](RightHalf)),
    print hex(ExpansionPerSbox[4](RightHalf)),
    print hex(ExpansionPerSbox[5](RightHalf)),
    print hex(ExpansionPerSbox[6](RightHalf)),
    print hex(ExpansionPerSbox[7](RightHalf))

    print "Inverse permutation per S-box:",
    print hex(InversePermutationPerSbox[0](RightHalf)),
    print hex(InversePermutationPerSbox[1](RightHalf)),
    print hex(InversePermutationPerSbox[2](RightHalf)),
    print hex(InversePermutationPerSbox[3](RightHalf)),
    print hex(InversePermutationPerSbox[4](RightHalf)),
    print hex(InversePermutationPerSbox[5](RightHalf)),
    print hex(InversePermutationPerSbox[6](RightHalf)),
    print hex(InversePermutationPerSbox[7](RightHalf))


##############################################################################
# Entrypoint for self-testing

if __name__ == "__main__":
    testDesUtilities()
    #generateInversePShifts()
