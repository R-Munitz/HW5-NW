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

    # assert that alignment matrix is correct
    expected_alignment_mat = np.array([[0., -10., -11., -12.],
                                        [-10.,   5.,  -5., -6.],
                                        [-11.,  -5.,   4.,  -6.],
                                        [-12., -6.,   0.,   5.],
                                        [-13.,  -7.,  -5.,   5.]])
                                            

    alignment_mat = nw._align_matrix
    assert np.array_equal(alignment_mat, expected_alignment_mat)

    #assert that gap matrices are correct
    gapA_mat = nw._gapA_matrix
    gapB_mat = nw._gapB_matrix

    # assert that gapA matrix is correct
    expected_gapA_mat = np.array([[-np.inf, -np.inf, -np.inf, -np.inf],
                                    [-10., -20., -21., -22.],
                                    [-11.,  -5., -15., -16.],
                                    [-12.,  -6.,  -6., -16.],
                                    [-13.,  -7.,  -7.,  -5.]])

    assert np.array_equal(gapA_mat, expected_gapA_mat)

    # assert that gapB matrix is correct
    expected_gapB_mat = np.array([[-np.inf, -10., -11., -12.],
                                    [-np.inf, -20.,  -5.,  -6.],
                                    [-np.inf, -21., -15.,  -6.],
                                    [-np.inf, -22., -16., -10.],
                                    [-np.inf, -23., -17., -15.]])

    assert np.array_equal(gapB_mat, expected_gapB_mat)

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

    # assert that backtrace matrix is correct
    backtrace_mat = nw._back
    # testing
    expected_backtrace_mat = np.array([[0., 0., 0., 0., 0., 0., 0., 0.],
                                        [0., 0., 2., 2., 2., 2., 2., 2.],
                                        [0., 1., 0., 0., 0., 0., 2., 0.],
                                        [0., 1., 1., 0., 0., 2., 2., 2.],
                                        [0., 1., 0., 1., 0., 0., 0., 2.],
                                        [0., 1., 0., 1., 0., 0., 0., 0.],
                                        [0., 1., 1., 0., 0., 1., 0., 0.],
                                        [0., 1., 1., 0., 0., 2., 2., 0.],
                                        [0., 1., 0., 1., 1., 0., 2., 2.],
                                        [0., 1., 0., 1., 1., 0., 0., 2.],
                                        [0., 1., 0., 1., 1., 1., 1., 0.]])
    assert np.array_equal(backtrace_mat, expected_backtrace_mat)

    # assert alignment score is correct
    assert score == 17

    # assert alignment seqs are correct
    assert align3 == "MAVHQLIRRP"
    assert align4 == "M---QLIRHP"

    pass
