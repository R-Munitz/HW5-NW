# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        #step 1- initialize matrices

        #matrix dimensions - add extra row and column for gaps
        row = len(seqB) + 1 #i seq B is vertical 
        col = len(seqA) + 1 #j seq A is horizontal
        
        #initialize matrices with zeros
        self._align_matrix = np.zeros((row, col))
        self._gapA_matrix = np.zeros((row, col))
        self._gapB_matrix = np.zeros((row, col))

        #initialize backtrace matrices with zeros
        self._back = np.zeros((row, col))
        self._back_A = np.zeros((row, col))
        self._back_B = np.zeros((row, col))

        """
        note: will one-hot encode/dict to store up,left,diagonal directions
        0 - diagonal
        1 - up
        2 - left
        """

        #step 2 - fill out first row and column with gap penalties

        #fill out first row with gap penalties - aligning seq A with gaps
        #gap opening penalty, then gap extend penalty for every subsequent cell

        for j in range(1, col): #move right across cols in first row
            self._align_matrix[0][j] = self.gap_open + j * self.gap_extend
            self._gapA_matrix[0][j] = -np.inf #can't align seqA with gap in seqA
            self._back[0][j] = 2 #left move
        
        #fill out first column with gap penalties - aligning seq B with gaps
        #gap opening penalty, then gap extend penalty for every subsequent cell
        for i in range(1, row): #move down across rows in first column
            self._align_matrix[i][0] = self.gap_open + i * self.gap_extend
            self._gapB_matrix[i][0] = -np.inf #can't align seqB with gap in seqB
            self._back[i][0] = 1 #up move
        
        #step 3 - fill out matrix with alignment scores
        """
        notes: NW algorithm
        diagonal move - from (i-1,j-1) aka align seqA[j-1] with seqB [i-1] , match/mismatch. score = previous score + sub_matrix[seqA[j-1], seqB[i-1]]
        up move - from (i-1, j), align seqA[i-1] with a gap in seqB, add gap to seqB. score = previous score + gap penalty
        left move - from (i, j-1),align seqB[j-1] with a gap in seqA, add gap to seqA. score = previous score + gap penalty
        Take max of the three scores
        """
        #calculate alignment scores for each cell in the matrix
        for i in range(1, row): #skip first row - gap penalties already filled in
            for j in range(1, col): #skip first col - gap penalties already filled in
                #calculate score for each possible move

                #diagonal move - match or mismatch score
                diagonal_score = self._align_matrix[i-1][j-1] + self.sub_dict[(seqA[j-1], seqB[i-1])] #look up score in substitution matrix
                
          
                #up move - adding gap to seqB
                up_score = max(self._align_matrix[i-1][j] + self.gap_open, self._gapB_matrix[i-1][j] + self.gap_extend)
                #takes the max of two possible scores
                #if a gap is already open, then previous score + gap extend penalty > current gap penalty + extend gap penalty
                #if a gap was not yet open, then previous score + gap open penalty > current gap penalty (-inf) + extend gap penalty
                #take the max of the two scores = checking if gap is already open or not

                #left move - adding gap to seqA
                left_score = max(self._align_matrix[i][j-1]+ self.gap_open, self._gapA_matrix[i][j-1] + self.gap_extend)
                #same logic as above

                #take max of the three scores
                score = max(diagonal_score, up_score, left_score)

                #update alignment matrix with score
                self._align_matrix[i][j] = score
                

                #update backtrace matrix with direction taken to get max score - (optimize with dictionary?)
                #update gap matrices with gap penalties
                if score == diagonal_score:
                    self._back[i][j] = 0 #diagonal move
                elif score == up_score:
                    self._back[i][j] = 1   #up move
                    #add gap to seqB - store penalty in gap matrix 
                    #open or extend gap
                    if self._gapB_matrix[i-1][j] == -np.inf: # no gap yet
                        self._gapB_matrix[i][j] = self.gap_open #add open gap penalty
                    else: #gap open, add extend gap penalty
                        self._gapB_matrix[i][j] = self.gap_extend
                else:
                    self._back[i][j] = 2   #left move
                    #add gap to seqA - store penalty in gap matrix
                    #open or extend gap
                    if self._gapA_matrix[i][j-1] == -np.inf: # no gap yet
                        self._gapA_matrix[i][j] = self.gap_open #add open gap penalty
                    else: #gap open, add extend gap penalty
                        self._gapA_matrix[i][j] = self.gap_extend
  
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """

        #start at bottom right corner of matrix
        i = len(self._seqB)
        j = len(self._seqA)
        
        #initialize aligned sequences
        seqA_align_list = []
        seqB_align_list = []

        #alignment score
        self.alignment_score = self._align_matrix[i][j] 

        #backtrace through matrix
        while i > 0 or j > 0:
            #which move was taken to get to current cell?
            move = self._back[i][j]
            if move == 0: #diagonal move
                seqA_align_list.append(self._seqA[j-1]) #add last res of seqA to alignment
                seqB_align_list.append(self._seqB[i-1]) #add last res of seqB to alignment

                #move diagonally
                i -= 1
                j -= 1

            elif move == 1: #up move - #add gap to seqB
                seqA_align_list.append(self._seqA[j-1]) #add last res of seqA to alignment
                seqB_align_list.append("-") #add gap to seqB

                #move up
                i -= 1

            else: #left move - add gap to seqA
                seqA_align_list.append("-") #add gap to seqA
                seqB_align_list.append(self._seqB[i-1]) #add last res of seqB to alignment

                #move left
                j -= 1

        #reverse lists
        seqA_align_list.reverse()
        seqB_align_list.reverse()

        #convert list to string - needed?
        self.seqA_align = "".join(seqA_align_list)
        self.seqB_align = "".join(seqB_align_list)

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
