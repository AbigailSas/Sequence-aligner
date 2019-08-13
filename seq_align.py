import argparse
import numpy as np
from itertools import groupby
import align_global
import local
import overlap

NAME = 0
SEQ = 1
TOP = 0
BOTTOM = 1

def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()

    parsed_a = fastaread(command_args.seq_a).__next__()
    a_name, a_seq = parsed_a[NAME], parsed_a[SEQ]
    parsed_b = fastaread(command_args.seq_b).__next__()
    b_name, b_seq = parsed_b[NAME], parsed_b[SEQ]
    score_matrix = np.genfromtxt(fname=command_args.score, delimiter="\t", skip_header=1, filling_values=1)[:, 1:]

    score = 0
    alignment = []
    if command_args.align_type == 'global':
        final_align = align_global.GlobalAlign(score_matrix, a_seq, b_seq)
        alignment = final_align.get_alignment()
        score = final_align.get_align_score()
    elif command_args.align_type == 'local':
        final_align = local.Local(score_matrix, a_seq, b_seq)
        alignment = final_align.get_alignment()
        score = final_align.get_align_score()
    elif command_args.align_type == 'overlap':
        final_align = overlap.Overlap(score_matrix, a_seq, b_seq)
        alignment = final_align.get_alignment()
        score = final_align.get_align_score()


    # print the best alignment and score

    split_size = 50
    length = len(alignment[TOP])
    top = [alignment[TOP][i:i+split_size] for i in range(0, length, split_size)]
    bottom = [alignment[BOTTOM][i:i+split_size] for i in range(0, length, split_size)]
    for i in range(len(top)):
        print(top[i] + "\n" + bottom[i] + "\n")
    print(command_args.align_type + ": " + str(score))

if __name__ == '__main__':
    main()
