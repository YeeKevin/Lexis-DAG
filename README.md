# COP5725 Advanced Databases: Lexis-DAG
Chu-Min Huang and Kevin Yee

Spring 2021

## Implementation of Lexis-DAG
Following the paper, "Lexis: An Optimization Framework for Discovering the Hierarchical Structure of Sequential Data"

## Project Tasks
- [x] Acquire Datasets
	+ Use the same or similar data as in original paper
	+ **Kevin Yee**
- [] Plan out design of the program
	+ A general blueprint for required functions and sections
- [] Implement features
	+ Create any necessary structures
	+ Greedy Algorithm
- [] Run Experiments
	+ Check results against original data

## Set-Up
Requries Bio-Python library
https://biopython.org/wiki/Download

GraphViz for printing of graphs
https://graphviz.org/download/
Use command to print pdf of DAG:
dot -Tpdf printedGraph.txt -o printedGraph.pdf

## Datasets

uniprot-proteome UP000002311.fasta
- Sequence of amino acids
- Redundancy reduced by UniProt's Proteome Redundancy Detector
- https://www.uniprot.org/uniprot/?query=proteome:UP000002311
- As of January 2021
- 6049 sequences

UP000002311_559292.fasta
- Sequence of amino acids
- Same as previous set except that for each gene, a representative seqence is selected
- 6049 sequences

bakers-yeast-6130.fasta
- From NCBI with search "Baker's Yeast" AND RefSeq filter enabled
- Sequence of nucleotides
- 6130 sequences
