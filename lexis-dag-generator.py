# Lexis-DAG Generator
# This program parses a fasta file and produces a Lexis-DAG for each target

from Bio import SeqIO               # processing .fasta format
from itertools import combinations  # enumerating substrings

# Global constants
MIN_OUT_DEGREE = 2      
MIN_SUBSTRING_LENGTH = 2
EMPTY_STRING = ""
SAMPLE_SIZE = 1

# Source amino acids
A = 'A'                 # alanine
R = 'R'                 # arginine
N = 'N'                 # asparagine
D = 'D'                 # aspartic acid
C = 'C'                 # cysteine
Q = 'Q'                 # glutamine
E = 'E'                 # glutamix acid
G = 'G'                 # glycine
H = 'H'                 # histidine
I = 'I'                 # isoleucine
L = 'L'                 # leucine
K = 'K'                 # lysine
M = 'M'                 # methionine
F = 'F'                 # phenylalanine
P = 'P'                 # proline
S = 'S'                 # serine
T = 'T'                 # threonine
W = 'W'                 # tryptophan
Y = 'Y'                 # tyrosine
V = 'V'                 # valine
B = 'B'                 # special case: asparagine/aspartic acid
Z = 'Z'                 # special case: glutamine/glutamic acid

acidArray = [A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z]

# class to represent DAG objects
class LexisDag:

    def __init(self):
        self.__target = EMPTY_STRING
        self.__source = dict()
        self.__intermediateNodes = dict()
        self.__edgeCost = 0

    # setters
    def setTarget(self, target):
        self.__target = target
    def setSource(self, source):
        self.__source = source
    def setEdgeCost(self, edgeCost):
        self.__edgeCost = edgeCost

    # getters
    def getTarget(self):
        return self.__target
    def getSource(self):
        return self.__source
    def getEdgeCost(self):
        return self.__edgeCost
        
    def __str__(self):
        dag_string = "\nTarget: \n" + str(self.__target) + \
                     "\nSource: \n" + str(self.__source) + \
                     "\nEdge Cost: " + str(self.__edgeCost)
        return dag_string

def main():

    # list of DAGs
    lexisDags = []
    iterator = 0

    #for seq_record in SeqIO.parse('dataset/UP000002311_559292.fasta', "fasta"):
    for seq_record in SeqIO.parse('dataset\\UP000002311_559292.fasta', "fasta"):
        # print(seq_record.id)

        # store the first 1500 sequences
        target = (seq_record.seq)
        initializeNewDag(target, lexisDags)
        iterator = iterator + 1
        if iterator == SAMPLE_SIZE:
            break

    # work with first target
    print(lexisDags[0])
    # all substrings
    target = str(lexisDags[0].getTarget())

    allSubstrings = repeatedSubstrings(target)
    print("Substrings: " + str(allSubstrings))
    print(len(str(allSubstrings)))

# create new Lexis-Dag and append to dagList
def initializeNewDag(target, dagList):
    newSequence = LexisDag()
    dagList.append(newSequence)
    # store target
    dagList[len(dagList) - 1].setTarget(target)
    # initial edge cost is length of target
    dagList[len(dagList) - 1].setEdgeCost(len(dagList[len(dagList) - 1].
        getTarget()))
    # current target
    target = dagList[len(dagList) - 1].getTarget()

    # set source nodes
    value = [pos for pos, char in enumerate(target) if char == R]
    dagList[len(dagList) - 1].setSource(
                  [{A : [pos for pos, char in enumerate(target) if char == A]},
                   {R : [pos for pos, char in enumerate(target) if char == R]},
                   {N : [pos for pos, char in enumerate(target) if char == N]},
                   {D : [pos for pos, char in enumerate(target) if char == D]},
                   {C : [pos for pos, char in enumerate(target) if char == C]},
                   {Q : [pos for pos, char in enumerate(target) if char == Q]},
                   {E : [pos for pos, char in enumerate(target) if char == E]},
                   {G : [pos for pos, char in enumerate(target) if char == G]},
                   {H : [pos for pos, char in enumerate(target) if char == H]},
                   {I : [pos for pos, char in enumerate(target) if char == I]},
                   {L : [pos for pos, char in enumerate(target) if char == L]},
                   {K : [pos for pos, char in enumerate(target) if char == K]},
                   {M : [pos for pos, char in enumerate(target) if char == M]},
                   {F : [pos for pos, char in enumerate(target) if char == F]},
                   {P : [pos for pos, char in enumerate(target) if char == P]},
                   {S : [pos for pos, char in enumerate(target) if char == S]},
                   {T : [pos for pos, char in enumerate(target) if char == T]},
                   {W : [pos for pos, char in enumerate(target) if char == W]},
                   {Y : [pos for pos, char in enumerate(target) if char == Y]},
                   {V : [pos for pos, char in enumerate(target) if char == V]},
                   {B : [pos for pos, char in enumerate(target) if char == B]},
                   {Z : [pos for pos, char in enumerate(target) if char == Z]},
                   ])

# return list of repeated substrings of length >= 2
def repeatedSubstrings(target):
    allSubstrings = [target[x:y] for x, y in combinations(
        range(len(target)+1), r = 2)]

    #TODO remove substrings of length 1

    #TODO remove non-repeated substrings

    return str(allSubstrings)

main()