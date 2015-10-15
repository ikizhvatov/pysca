'''
DES transformations required for DPA with conditional averaging

Uses minor chunks of code from pyDES-2.0.1 and DPA contest v1 DES example.

TODO: rewrite in Cython or in C using cyclic shifts and other natural bitwise
      operations; look at DES implementation in libtomcrypt as an example.

Started by Ilya on 2014-11-25
'''

from operator import sub


##############################################################################
# Core functionality

''' Bit permutations '''
def permuteBits(x, permutation, inputLength):
    ''' Permutes bits of x given a permutation table. Assumes that permutation table is 0-offset.
        The input bitlength is a parameter
        The output bitlength is determined by the permutation table
    '''
    result = 0L
    for i in range(0, len(permutation)):
        result = ((result << 1) |
                  ((x >> (inputLength - 1 - permutation[i])) & 1))
    return result

# These are bit permutaions to be used with permuteBits above,
# not lookup tables.
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
RoundPermutation = [
    15,  6, 19, 20, 28, 11, 27, 16,
    0,  14, 22, 25,  4, 17, 30,  9,
    1,   7, 23, 13, 31, 26,  2,  8,
    18, 12, 29,  5, 21, 10,  3, 24
    ]
PC1Permutation = [
    56, 48, 40, 32, 24, 16,  8,
     0, 57, 49, 41, 33, 25, 17,
     9,  1, 58, 50, 42, 34, 26,
    18, 10,  2, 59, 51, 43, 35,
    62, 54, 46, 38, 30, 22, 14,
     6, 61, 53, 45, 37, 29, 21,
    13,  5, 60, 52, 44, 36, 28,
    20, 12,  4, 27, 19, 11,  3
    ]
PC2Permutation = [
    13, 16, 10, 23,  0,  4,
     2, 27, 14,  5, 20,  9,
    22, 18, 11,  3, 25,  7,
    15,  6, 26, 19, 12,  1,
    40, 51, 30, 36, 46, 54,
    29, 39, 50, 44, 32, 47,
    43, 48, 38, 55, 33, 52,
    45, 41, 49, 35, 28, 31
    ]

''' S-box '''
def sBox(m, x):
    row = ((x & 0x20) >> 4) ^ (x & 1) 
    col = (x & 0x1e) >> 1
    return SBoxLUT[m][16 * row + col]

# This is a 6-to-4 bits lookup table. It is not directly
# addressable with S-box input but requires a row-col transform,
# see sBox() above.
SBoxLUT = [
    # S1
    [14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7,
     0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8,
     4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0,
     15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13],

    # S2
    [15, 1, 8, 14, 6, 11, 3, 4, 9, 7, 2, 13, 12, 0, 5, 10,
     3, 13, 4, 7, 15, 2, 8, 14, 12, 0, 1, 10, 6, 9, 11, 5,
     0, 14, 7, 11, 10, 4, 13, 1, 5, 8, 12, 6, 9, 3, 2, 15,
     13, 8, 10, 1, 3, 15, 4, 2, 11, 6, 7, 12, 0, 5, 14, 9],

    # S3
    [10, 0, 9, 14, 6, 3, 15, 5, 1, 13, 12, 7, 11, 4, 2, 8,
     13, 7, 0, 9, 3, 4, 6, 10, 2, 8, 5, 14, 12, 11, 15, 1,
     13, 6, 4, 9, 8, 15, 3, 0, 11, 1, 2, 12, 5, 10, 14, 7,
     1, 10, 13, 0, 6, 9, 8, 7, 4, 15, 14, 3, 11, 5, 2, 12],

    # S4
    [7, 13, 14, 3, 0, 6, 9, 10, 1, 2, 8, 5, 11, 12, 4, 15,
     13, 8, 11, 5, 6, 15, 0, 3, 4, 7, 2, 12, 1, 10, 14, 9,
     10, 6, 9, 0, 12, 11, 7, 13, 15, 1, 3, 14, 5, 2, 8, 4,
     3, 15, 0, 6, 10, 1, 13, 8, 9, 4, 5, 11, 12, 7, 2, 14],

    # S5
    [2, 12, 4, 1, 7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9,
     14, 11, 2, 12, 4, 7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6,
     4, 2, 1, 11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14,
     11, 8, 12, 7, 1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3],

    # S6
    [12, 1, 10, 15, 9, 2, 6, 8, 0, 13, 3, 4, 14, 7, 5, 11,
     10, 15, 4, 2, 7, 12, 9, 5, 6, 1, 13, 14, 0, 11, 3, 8,
     9, 14, 15, 5, 2, 8, 12, 3, 7, 0, 4, 10, 1, 13, 11, 6,
     4, 3, 2, 12, 9, 5, 15, 10, 11, 14, 1, 7, 6, 0, 8, 13],

    # S7
    [4, 11, 2, 14, 15, 0, 8, 13, 3, 12, 9, 7, 5, 10, 6, 1,
     13, 0, 11, 7, 4, 9, 1, 10, 14, 3, 5, 12, 2, 15, 8, 6,
     1, 4, 11, 13, 12, 3, 7, 14, 10, 15, 6, 8, 0, 5, 9, 2,
     6, 11, 13, 8, 1, 4, 10, 7, 9, 5, 0, 15, 14, 2, 3, 12],

    # S8
    [13, 2, 8, 4, 6, 15, 11, 1, 10, 9, 3, 14, 5, 0, 12, 7,
     1, 15, 13, 8, 10, 3, 7, 4, 12, 5, 6, 11, 0, 14, 9, 2,
     7, 11, 4, 1, 9, 12, 14, 2, 0, 6, 10, 13, 15, 3, 5, 8,
     2, 1, 14, 7, 4, 10, 8, 13, 15, 12, 9, 0, 3, 5, 6, 11]
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

''' Key expansion '''
def computeRoundKeys(key, numberOfRounds):

    keyShifts = [1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1]
    mask28 = 0xfffffff

    # rotate left modulo 28 bits
    rol28 = lambda x, n: ((x << n) & mask28) | ((x & mask28) >> (28 - n))

    permutedKey = permuteBits(key, PC1Permutation, 64) 
    
    l = (permutedKey >> 28) & mask28
    r = permutedKey & mask28

    roundKeys = []
    for i in range(numberOfRounds):
        l = rol28(l, keyShifts[i])
        r = rol28(r, keyShifts[i])
        lr = (l << 28) ^ r
        roundKey = permuteBits(lr, PC2Permutation, 56)
        roundKeys.append(roundKey)
        
    return roundKeys

''' Return n-th 6-bit chunk of the 48-bit round key '''
def roundKeyChunk(roundKey, n):
    return (roundKey >> (42 - 6 * n)) & 0x3f


##############################################################################
# Tandem of functions for round in xor out intermediate. First functions
# obtains value for conditional averaging. Second function obtains the value
# of the target variable from that

def roundXOR_valueForAveraging(input, sBoxNumber):
    ''' Compute the value for conditional averaging from input, for a given
        S-box number '''

    # prepare the first round input halves
    permutedInput = permuteBits(input, InitialPermutation, 64)
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
    SBoxOut = sBox(sBoxNumber, SBoxIn)
    RoundInXorOutPerSBox = SBoxOut ^ y

    return RoundInXorOutPerSBox

# Both merged into none, for attack without conditional averaging
def roundXOR_allInOne(input, keyChunk, sBoxNumber):

    # prepare the first round input halves
    permutedInput = permuteBits(input, InitialPermutation, 64)
    rightHalf = permutedInput & 0xFFFFFFFF
    leftHalf = permutedInput >> 32

    # get S-box output
    a = ExpansionPerSbox[sBoxNumber](rightHalf) # returns 6 bits of S-box input
    SBoxIn = a ^ keyChunk
    SBoxOut = sBox(sBoxNumber, SBoxIn)
    
    # gather the input bits that need to be XORed with the S-box output
    b = InversePermutationPerSbox[sBoxNumber](rightHalf ^ leftHalf)
    
    # compute the XOR
    RoundInXorOutPerSBox = SBoxOut ^ b
    
    return RoundInXorOutPerSBox


##############################################################################
# Self-creators

def generateInversePermutationPerSbox():
    ''' Helper used to generate the shifts. In the output, negative values should be manually replaced by a left shift! '''
    print '\n--- generateInversePermutationPerSbox ---'

    initialPositionsPerSbox = [
        [ 8, 16, 22, 30],
        [12, 27,  1, 17],
        [23, 15, 29,  5],
        [25, 19,  9,  0],
        [ 7, 13, 24,  2],
        [ 3, 28, 10, 18],
        [31, 11, 21,  6],
        [ 4, 26, 14, 20]
        ]
    finalPositions = [28, 29, 30, 31]

    for group in initialPositionsPerSbox:
        shifts = map(sub, finalPositions, group) # element-wise list subtraction
        print "((x >> %d) & 8) | ((x >> %d) & 4) | ((x >> %d) & 2) | ((x >> %d) & 1)" % (shifts[0], shifts[1], shifts[2], shifts[3])


##############################################################################
# Self-tests

def testDesUtilities():
    ''' Dump the state of the first round to compare against a reference implementation.
        Compare the inverse round permutation against the forward one.
        The output should look like:

        --- testDesUtilites ---
        L  : 0x59e0bc92L
        R  : 0xa69230c8L
        RK0: 0x8805bc20c812L
        Rt : 0x50d4a41a1651L
        Rtk: 0xd8d1183ade43L
        z  : 0x789b6fef
        zp : 0x9c7eafebL
        Testing the inverse permutation
        z' : 0x789b6fefL
        Success!
    '''
    print '\n--- testDesUtilites ---'

    # data from the first trace in TC8 PA training traceset
    key        = 0x8a7400a03230da28L
    plaintext  = 0x40a184466d9c52b7L
    ciphertext = 0x1cb5ca37b8a7a388L

    # key schedule
    roundKeys = computeRoundKeys(key, 16)
    k = roundKeys[0]

    # prepare the first round input halves (checked)
    permutedInput = permuteBits(plaintext, InitialPermutation, 64)
    rightHalf = permutedInput & 0xFFFFFFFF
    leftHalf = permutedInput >> 32
    print 'L  : ' + hex(leftHalf)
    print 'R  : ' + hex(rightHalf)
    print 'RK0: ' + hex(k)

    #  expansion (checked)
    Rt = 0L
    for i in range(0, 8):
       a = ExpansionPerSbox[i](rightHalf)
       Rt = (Rt << 6) ^ a;
    print 'Rt : ' + hex(Rt)

    # key addition
    Rt = Rt ^ k
    print 'Rtk: ' + hex(Rt)

    # S-boxes
    z = 0L
    for i in range(0, 8):
        z ^= (sBox(7 - i, Rt & 0x3f) << (i * 4))
        Rt = Rt >> 6
    print 'z  : ' + hex(z)

    # permutation
    zp = permuteBits(z, RoundPermutation, 32)
    print 'zp : ' + hex(zp)

    # testing the inverse permutation
    print 'Testing the inverse permutation'
    zb = 0L
    for i in range(0, 8):
        zb ^= InversePermutationPerSbox[i](zp) << ((7 - i) * 4)
    print "z' : " + hex(zb)
    if (zb == z):
        print 'Success!'
    else:
        print 'Fail!'


def dumpRoundKeys():
    ''' Dump all round keys '''
    print '\n--- dumpRoundKeys ---'

    key = 0x8a7400a03230da28L
    roundKeys = computeRoundKeys(key, 16)
    print 'Key  : ' + format(key, '#018x')
    for i in range(16):
        print 'RK' + format(i, '02d') + ' : ' + format(roundKeys[i], '#014x'),
        print '[',
        for j in range(8):
            print format(roundKeyChunk(roundKeys[i], j), '#04x'),
        print ']'

def dumpMiscValues():
    ''' Print out the values, just in case '''
    print '\n--- dumpMiscValues ---'

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
    print hex(permuteBits(Input, InitialPermutation, 64))

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
    dumpRoundKeys()
    dumpMiscValues()
    generateInversePermutationPerSbox()
