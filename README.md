# Sequence-aligner
A dna sequence aligner that give 2 sequences, a choice of alignment option and score function prints optimal alignment. 
The algorithms used are dynamic programming algorithms. 

Input: 
1. 2 fasta files with 1 sequence each, over the alphabet A,C,T,G
2. Scoring matrix - score for match, deletion or mismatch, given as a tsv file. 

Output: 
Printing the alignment of 2 given sequences X,Y,  according to the score function. User needs to choose one of the following alignment options: 
1. global alignment - the algorithm searches for the best match between X and Y such that all bases xi are aligned to some yj or to gaps. 
2. local alignment - the algorithm searches for the best alignment of substrings of X and Y. 
3. overlap alignment - this is usefull if one sequence is a substring of the other, or if they overlap (the end of one sequence with the begening of the other). the algorithm searches for the best global aligmemt with no penelty for gaps in the beggening or the end.

The program can be run by writing the command: 
"python3 seq_align.py a.fasta b.fasta --align_type global --score score_matrix.tsv"

