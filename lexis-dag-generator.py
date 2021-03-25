from Bio import SeqIO

#for seq_record in SeqIO.parse('dataset/UP000002311_559292.fasta', "fasta"):
# windows path
for seq_record in SeqIO.parse('dataset\\UP000002311_559292.fasta', "fasta"):
	# print(seq_record.id)

	# access each amino acid sequence
	print(seq_record.seq)