# Lexis-DAG Generator
# This program parses a fasta file and produces a Lexis-DAG for each target

from Bio import SeqIO

# Global constants
MIN_OUT_DEGREE = 2		
MIN_SUBSTRING_LENGTH = 2
EMPTY_STRING = ""
SAMPLE_SIZE = 1500

# Source amino acids
A = 'A'					# alanine
R = 'R'					# arginine
N = 'N'					# asparagine
D = 'D'					# aspartic acid
C = 'C'					# cysteine
Q = 'Q'					# glutamine
E = 'E'					# glutamix acid
G = 'G'					# glycine
H = 'H'					# histidine
I = 'I'					# isoleucine
L = 'L'					# leucine
K = 'K'					# lysine
M = 'M'					# methionine
F = 'F'					# phenylalanine
P = 'P'					# proline
S = 'S'					# serine
T = 'T'					# threonine
W = 'W'					# tryptophan
Y = 'Y'					# tyrosine
V = 'V'					# valine
B = 'B'					# special case: asparagine/aspartic acid
Z = 'Z'					# special case: glutamine/glutamic acid

# class to represent DAG objects
class LexisDag:

	def __init(self):
		self.__target = EMPTY_STRING
		self.__source = []
		self.__intermediateNodes = []
		self.__edgeCost = 0

	# setters
	def setTarget(self, target):
		self.__target = target
	def setEdgeCost(self, edgeCost):
		self.__edgeCost = edgeCost

	# getters
	def getTarget(self):
		return self.__target
	def getEdgeCost(self):
		return self.__edgeCost

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
	print(lexisDags[0].getTarget())
	print(lexisDags[0].getEdgeCost())


def initializeNewDag(target, dagList):
	newSequence = LexisDag()
	dagList.append(newSequence)
	# store target
	dagList[len(dagList) - 1].setTarget(target)
	# initial edge cost is length of target
	dagList[len(dagList) - 1].setEdgeCost(len(dagList[len(dagList) - 1].getTarget()))


	# TODO add in source


main()