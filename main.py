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
    
    seq_species_dict = {
        gg_seq:"Gallus gallus",
        mm_seq:"Mus musculus",
        br_seq:"Balaeniceps rex",
        tt_seq:"Tursiops truncatus"
    }
    
    """
    test_seq = "AGT"
    test_seq_b = "ALGT"
    
    score, hs_seq_align, seq_B_align = nm.align(test_seq, test_seq_b)
    print("Alignment score: ", score)
    print("hs_seq_align: ", hs_seq_align)
    print("seq_B_align: ", seq_B_align)
    """
    
    for seq in sequences_to_align:
        #align seq to hs_seq
        score, hs_seq_align, seq_b_align = nm.align(hs_seq, seq)
        seq_score_dict[seq] = score


    #print species in order of most similar to human BRD
    sorted_seq_score_dict = dict(sorted(seq_score_dict.items(), key=lambda item: item[1], reverse=True))
    for key in sorted_seq_score_dict.keys():
        print(f"Species: ", seq_species_dict[key], "; Alignment score: ", sorted_seq_score_dict[key])


def troubleshooting():
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    nw = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    score, align1, align2 = nw.align(seq1, seq2)

    print("alignment matrix: ", nw._align_matrix)
    print("gapA matrix: ", nw._gapA_matrix)
    print("gapB matrix: ", nw._gapB_matrix)
    
     


    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    score, align3, align4 = nw.align(seq3, seq4)

    print("score: ", score)
    print("align3: ", align3)
    print("align4: ", align4)

    print("nw backtrace: ", nw._back)


if __name__ == "__main__":
    main()
    troubleshooting()
