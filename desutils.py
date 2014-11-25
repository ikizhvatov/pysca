'''
DES transformations required for DPA with conditional averaging

TODO: rewrite in Cython or in C using cyclic shifts and other natural bitwise
      operations; look at DES implementation in libtomcrypt as an example.

Started by Ilya on 2014-11-25
'''

import numpy as np


##############################################################################
# Core functionality

def permuteBits(x, permutation):
    ''' Permutes bits of x given a permutation table. Implemented as in Inspector 4 code. Assumes that permutation table is 0-offset. '''

    result = 0L
    for i in range(0, len(permutation)):
        result = (result << 1) | ((x >> (len(permutation) - 1 - permutation[i])) & 1)
    return result

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


''' Extension permutation, as a list of function per S-box.
    x is supposed to be a 32-bit wide integer. If not, function will work incorrectly.
    Calling example: ExtensionPermutationPerSbox[4](x) - retrive the part of E(x) corresponding to S-box 4'''
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

'''Inverse P permutation as a list of bit-gathering functions per S-box. Generated using a helper below.'''
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


def valueForAveraging_RoundXOR(input, sBoxNumber):
    ''' Compute the value for conditional averaging from input, for a given S-box number '''

    # prepare the first round input halves
    permutedInput = permuteBits(input, InitialPermutation)
    rightHalf = permutedInput & 0xFFFFFFFF
    leftHalf = permutedInput >> 32

    # get first part: the input chunks per S-box
    # concatenate; TODO: consider returning as some structure, if can average based on the structre value
    # get second part: the bits from XOR of left and right input 
    a = ExpansionPerSbox[sBoxNumber](rightHalf)
    b = InversePermutationPerSbox[sBoxNumber](rightHalf ^ leftHalf)
    r = (a << 4) ^ b

    return r

# TODO compute the target variable using the key chunk and the value for conditional averaging


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

def test():
    ''' Unit test for DES utilities '''
    # TODO check against a test vector generated in Inspector 4 or elsewhere

    input = 0xA76DB873C63FE078
    print hex(input)

    print "The full thing: "
    print hex(valueForAveraging_RoundXOR(input, 0))
    print hex(valueForAveraging_RoundXOR(input, 1))
    print hex(valueForAveraging_RoundXOR(input, 2))
    print hex(valueForAveraging_RoundXOR(input, 3))
    print hex(valueForAveraging_RoundXOR(input, 4))
    print hex(valueForAveraging_RoundXOR(input, 5))
    print hex(valueForAveraging_RoundXOR(input, 6))
    print hex(valueForAveraging_RoundXOR(input, 7))


    print "Initial permutation"
    print hex(permuteBits(input, InitialPermutation))

    rightHalf = input & 0xFFFFFFFF

    print "Expansion per S-box:",
    print hex(ExpansionPerSbox[0](rightHalf)),
    print hex(ExpansionPerSbox[1](rightHalf)),
    print hex(ExpansionPerSbox[2](rightHalf)),
    print hex(ExpansionPerSbox[3](rightHalf)),
    print hex(ExpansionPerSbox[4](rightHalf)),
    print hex(ExpansionPerSbox[5](rightHalf)),
    print hex(ExpansionPerSbox[6](rightHalf)),
    print hex(ExpansionPerSbox[7](rightHalf))

    print "Inverse permutation per S-box:",
    print hex(InversePermutationPerSbox[0](rightHalf)),
    print hex(InversePermutationPerSbox[1](rightHalf)),
    print hex(InversePermutationPerSbox[2](rightHalf)),
    print hex(InversePermutationPerSbox[3](rightHalf)),
    print hex(InversePermutationPerSbox[4](rightHalf)),
    print hex(InversePermutationPerSbox[5](rightHalf)),
    print hex(InversePermutationPerSbox[6](rightHalf)),
    print hex(InversePermutationPerSbox[7](rightHalf))


##############################################################################
# Entrypoint for self-testing

if __name__ == "__main__":
    test()
    #generateInversePShifts()