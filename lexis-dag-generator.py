from Bio import SeqIO

# path is for windows environment
# dataset/bakers-yeast-6130.fasta
for seq_record in SeqIO.parse('dataset\\bakers-yeast-6130.fasta', "fasta"):
	print(seq_record.id)
	print(seq_record.seq)