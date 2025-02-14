# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    nw = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, align1, align2 = nw.align(seq1, seq2)

    #assert that alignment matrix is correct
    #expected_alignment_mat = [[  0., -10., -11., -12.], #double check this is correct?
                                #[-10.,   5.,  -5.,  -6.],
                                #[-11.,  -5.,   4.,   2.],
                                #[-12.,  -6.,  10.,   9.]]
    #testing
    expected_alignment_mat = [
    [  0., -10., -11., -12.],  
    [-10.,   5.,  -5.,  -6.],  
    [-11.,  -5.,   4.,   2.],  
    [-12.,  -6.,  10.,   9.],  
    [-13.,  -7.,   1.,  15.]  ]
                                
    alignment_mat = nw._align_matrix
    assert np.array_equal(alignment_mat, np.array(expected_alignment_mat))

    #note, expected matrices are just placeholders, need to check if they are correct
    gapA_mat = nw.gapA_matrix
    gapB_mat = nw.gapB_matrix

    #assert that gapA matrix is correct
    #testing
    expected_gapA_mat = [
    [  0., -10., -11., -12.],  
    [-10., -np.inf, -20., -21.],  
    [-11., -np.inf, -15., -16.],  
    [-12., -np.inf, -16., -17.],  
    [-13., -np.inf, -17., -18.]]
    
    assert np.array_equal(gapA_mat, np.array(expected_gapA_mat))

    #assert that gapB matrix is correct
    #testing
    expected_gapB_mat = [
    [  0., -10., -11., -12.],  
    [-10., -np.inf, -np.inf, -np.inf],  
    [-11., -20., -15., -16.],  
    [-12., -21., -16., -17.],  
    [-13., -22., -17., -18.] ]
    assert np.array_equal(gapB_mat, np.array(expected_gapB_mat))


    pass
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    score, align3, align4 = nw.align(seq3, seq4)

    #assert that backtrace matrix is correct
    backtrace_mat = nw._back
    #testing
    expected_backtrace_mat = [
    [0, 2, 2, 2, 2, 2, 2, 2],
    [1, 0, 2, 2, 2, 2, 2, 2],
    [1, 1, 2, 2, 2, 2, 2, 2],
    [1, 1, 1, 2, 2, 2, 2, 2],
    [1, 1, 1, 1, 2, 2, 0, 2],
    [1, 1, 0, 1, 2, 2, 2, 2],
    [1, 1, 1, 0, 2, 2, 2, 2],
    [1, 1, 1, 1, 0, 2, 2, 2],
    [1, 1, 1, 1, 1, 0, 2, 2],
    [1, 1, 1, 1, 1, 0, 2, 2],
    [1, 1, 1, 1, 1, 1, 2, 0]
]
    assert np.array_equal(backtrace_mat, np.array(expected_backtrace_mat))

    #assert alignment score is correct
    assert score == 17

    #assert alignment seqs are correct
    assert align3 == "MAVHQLIRRP"
    assert align4 == "M---QLIRHP"

    pass




