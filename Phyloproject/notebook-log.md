Course project

My research aims to reconstruct a complete phylogeny of the neotropical clades of the orchid genus Vanilla (subgen. Vanilla and sect. Xanata). Since I don't have my own data yet, I'll use accessions posted in GenBank, the DNA regions chosen are matK and rbcL.


[Feb23]

I had several problems trying to install clustalw. Finally I used this command and it worked!

	conda create -n clustalw2 -c biobuilds -y clustalw

To activate clustalw the command is:

	conda activate clustalw2

To run the alignment first you need to be located inside the folder that has the fasta file, then run the following command:

	clustalw2 -ALIGN -INFILE=primatesAA.fasta -OUTFILE=primatesAA-aligned.fasta -OUTPUT=PHYLIP

To merge various .fasta files use "cat" command. In this case "> 1.fasta" means that all sequences will me merge on that file, this will overwrite whatever is inside that file

	cat 2.fasta 3.fasta 4.fasta 5.fasta 6.fasta 7.fasta 8.fasta 9.fasta 10.fasta 11.fasta 12.fasta 13.fasta 14.fasta 15.fasta 16.fasta 17.fasta 18.fasta 19.fasta 20.fasta 21.fasta >1.fasta

