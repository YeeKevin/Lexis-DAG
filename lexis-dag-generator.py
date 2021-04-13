# Lexis-DAG Generator
# This program parses a fasta file and produces a Lexis-DAG for each target

from Bio import SeqIO               # processing .fasta format
from itertools import combinations  # enumerating substrings
import re

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

# set a variable for the best costed result
bestSubList = []

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
    def setIntermediateNodes(self, intermediateNodes):
        self.__intermediateNodes = intermediateNodes
    def setEdgeCost(self, edgeCost):
        self.__edgeCost = edgeCost

    # getters
    def getTarget(self):
        return self.__target
    def getSource(self):
        return self.__source
    def getIntermediateNodes(self):
        return self.__intermediateNodes
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

    for seq_record in SeqIO.parse('dataset/UP000002311_559292.fasta', "fasta"):
    # for seq_record in SeqIO.parse('dataset\\UP000002311_559292.fasta', "fasta"):
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
    source = lexisDags[0].getSource()
    createPrintedGraph(target, source, "initialGraph.txt")

    # allSubstrings = repeatedSubstrings(target)
    # print("Substrings: " + str(allSubstrings))
    # print(len(str(allSubstrings)))
    
    # set test dataset as paper did, can put any dataset string in it
    getSubstrings(target, 0)
    print(bestSubList)

    # longest substring for experiment
    # print(longestSubstringHeuristic(bestSubList))

    # Lexis-G
    # parallel list of savedCost values
    savedCostList = calcSavedCost(bestSubList, target, [])
    print(savedCostList)

    # sorting from most saved cost
    zippedSubstrings = zip(savedCostList, bestSubList)
    sortedZippedSubstrings = sorted(zippedSubstrings, reverse=True)
    # print(sortedZippedSubstrings)

    # substrings sorted from highest to lowest savedCost
    sortedList = [element for _,element in sortedZippedSubstrings]
    print(sortedList)

    # TODO insert intermediate nodes
    insertIntermediateNodes(lexisDags[0], sortedList, target)

    print(lexisDags[0].getIntermediateNodes())
    print(lexisDags[0].getEdgeCost())

                
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
                   # {B : [pos for pos, char in enumerate(target) if char == B]},
                   # {Z : [pos for pos, char in enumerate(target) if char == Z]},
                   ])

# return list of repeated substrings of length >= 2
def repeatedSubstrings(target):
    allSubstrings = [target[x:y] for x, y in combinations(
        range(len(target)+1), r = 2)]

    #TODO remove substrings of length 1

    #TODO remove non-repeated substrings

    return str(allSubstrings)

# function that set best costed(1st layer) as an index in order to find best costed(2nd to last layer)
# paper target="aabcaabdaabc", best costed(1st layer) = "aab", continue to find if there is any best costed
# therefore, find best the 2nd best costed(2nd layer = "aabc"
def getSubstrings(target, index):
    Lt = len(target)
    # stemp = the best costed one
    stemp = ""
    ctemp = 0
    result = 0
    templist = []
    maxstep = int(Lt/2) + 1
    for step in range(2, maxstep):
        for i in range(0, Lt):
            if i+step == Lt+1:
               break
            sub = target[i: i+step]
            templist.append(sub)
            R = 0
            slen = 1
            j = 0
            while j < Lt:
                sub2 = target[j: j+step]
                if sub == sub2:
                    R+=1
                    slen = step
                else:
                    slen = 1
                j+=slen
            if R == 1:
               continue
            else:
                if (R-1)*step > ctemp:
                   stemp = sub
                   ctemp = (R-1)*step
    if  stemp == "":
        return bestSubList
    else:
        if index == 0:
            bestSubList.append(stemp)
        # if there is only one best costed substring, it will directly displayed
        else:
            t = stemp.replace(str(index-1), bestSubList[index-1])
            bestSubList.append(t)
        # "aab"(1st best costed) now is represented as index(0) in order to code easily, 
        # so if there is one more layer for best costed, combining them here, 
        # "0c" = "aab"(0) + "c". So, the final best costed one is "aabc".
        target = target.replace(stemp, str(index))
        index += 1
        getSubstrings(target, index)

# create file for printing of graph
# use command dot -Tpdf initialGraph.txt -o initialGraph.pdf
def createPrintedGraph(target, source, filename):

    f = open(filename, "w")
    # file prologue
    f.write("digraph G {\n")

    # iterate through source nodes
    i = 0
    for aminoAcid in source:
        # print(str(aminoAcid))
        currentKey = str(next(iter(aminoAcid)))
        f.write("    " + currentKey +  "->" + target + "[label=\"" + 
            str(aminoAcid[currentKey]) + "\"]"  + ";" + "\n")

    # iterate through intermediate nodes

    # file epilogue
    f.write("}")
    f.close()

# sort list from longest to shortest substrings
def longestSubstringHeuristic(unsortedSubstrings):
    sortedSubstrings = sorted(unsortedSubstrings, key=len, reverse=True)
    # sortedSubstrings.reverse()
    return sortedSubstrings

# calculate savedCost for a list of substrings
def calcSavedCost(bestSubList, target, intermediateNodes):
    savedCostList = []
    for substring in bestSubList:
        substringOccurences = target.count(substring)
        # print(substringOccurences)
        for intNode in intermediateNodes:
            if intNode.count(substring):
                substringOccurences = substringOccurences 
                + intNode.count(substring)
        # savedCost is occurences multiplied by length of substring - 1
        savedCost = substringOccurences * (len(substring) - 1)
        savedCostList.append(savedCost)
    # returns a parallel list of values
    return savedCostList

def insertIntermediateNodes(lexisDagObject, sortedList, target):

    # create node
    intList = []

    # pointing to target
    for node in sortedList:
        occurences = [m.start() for m in re.finditer(node, target)]
        entry = { node : occurences}
        intList.append(entry)
        # original edge cost of substring
        originalCost = len(occurences) * (len(node))
        # cost of adding the intermediate node
        newCost = len(occurences) + (len(node))
        difference = originalCost - newCost
        # update edge cost for DAG
        lexisDagObject.setEdgeCost(lexisDagObject.getEdgeCost() - difference)

        # edge cost in other int ndoes
        # otherIntNodes = lexisDagObject.getIntermediateNodes()


    # print(intList)
    
    # update edge cost (target)
    # add edges to construct node as edges to the target

    # udpate edge cost (current int nodes)


    # TODO edges to other int nodes
    lexisDagObject.setIntermediateNodes(intList)


    return


main()