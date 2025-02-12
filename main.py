# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch


def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    #align all provided species BRD2 sequences to the human BRD2 sequence and print the species in order of most similar to least similar with respect to human BRD2.
    #print the alignment scores corresponding to each species alignment to the human BRD2 sequence.     
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nm = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    #align each seq to human BRD
    sequences_to_align = [gg_seq,mm_seq,br_seq,tt_seq]
    #dictionary of sequences:alignment scores
    seq_score_dict = {
        gg_seq:None,
        mm_seq:None,
        br_seq:None,
        tt_seq:None

    }

    for seq in sequences_to_align:
        #print seq
        #align seq to hs_seq
        score, hs_seq_align, seq_B_align = nm.align(hs_seq, seq)
        seq_score_dict[seq] = score

    #print in order of highest score to lowest score
    sorted_seq_score_dict = dict(sorted(seq_score_dict.items(), key=lambda item: item[1], reverse=True))
    for key in sorted_seq_score_dict.keys():
        print(key)
        print("Alignment score: ", sorted_seq_score_dict[key])

    

    pass

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    pass
    

if __name__ == "__main__":
    main()
