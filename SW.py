from enum import IntEnum
import time

# Assigning the constant values for the traceback
class Trace(IntEnum):
    STOP = 0
    LEFT = 1 
    UP = 2
    DIAGONAL = 3
    
# Scores - to be given by the main function
class Score(IntEnum):
    MATCH = 2
    MISMATCH = -1
    GAP = -1

# Implementing the Smith Waterman local alignment
def smith_waterman(seq1, seq2):
    # Generating the empty matrices for storing scores and tracing
    row = len(seq1) + 1
    col = len(seq2) + 1
    score_matrix = [[0]*col for i in range(row+1)]  
    tracing_matrix = [[0]*col for i in range(row+1)]  
    
    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)
    
    # Calculating the scores for all cells in the score_matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = Score.MATCH if seq1[i - 1] == seq2[j - 1] else Score.MISMATCH
            diagonal_score = score_matrix[i - 1][ j - 1] + match_value
            
            # Calculating the vertical gap score
            vertical_score = score_matrix[i - 1][j] + Score.GAP
            
            # Calculating the horizontal gap score
            horizontal_score = score_matrix[i][j - 1] + Score.GAP
            
            # Taking the highest score 
            score_matrix[i][j] = max(0, diagonal_score, vertical_score, horizontal_score)
            
            # Tracking where the cell's value is coming from    
            if score_matrix[i][j] == 0: 
                tracing_matrix[i][j] = Trace.STOP
                
            elif score_matrix[i][j] == horizontal_score: 
                tracing_matrix[i][j] = Trace.LEFT
                
            elif score_matrix[i][j] == vertical_score: 
                tracing_matrix[i][j] = Trace.UP
                
            elif score_matrix[i][j] == diagonal_score: 
                tracing_matrix[i][j] = Trace.DIAGONAL 
                
            # Tracking the cell with the maximum score
            if score_matrix[i][j] >= max_score:
                max_index = (i,j)
                max_score = score_matrix[i][j]
    
    # Initialising the variables for tracing
    aligned_seq1 = ""
    aligned_seq2 = ""   
    current_aligned_seq1 = ""   
    current_aligned_seq2 = ""  
    (max_i, max_j) = max_index
    
    # Tracing and computing the pathway with the local alignment
    while tracing_matrix[max_i][max_j] != Trace.STOP:
        if tracing_matrix[max_i][max_j] == Trace.DIAGONAL:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1
            
        elif tracing_matrix[max_i][max_j] == Trace.UP:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'
            max_i = max_i - 1    
            
        elif tracing_matrix[max_i][max_j] == Trace.LEFT:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]
            max_j = max_j - 1
            
        aligned_seq1 = aligned_seq1 + current_aligned_seq1
        aligned_seq2 = aligned_seq2 + current_aligned_seq2
    
    # Reversing the order of the sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]
    
    return aligned_seq1, aligned_seq2

if __name__ == "__main__":
    t_start = time.time()

    R = "AAAGGCTGGGGACCACTGATCTAAATACACCAATAAAAAGAAAAAGATTGTAAGATTGGAGTTTAAAAGACCTGACTCTATACTGACCACAAATAAAAACC"
    Q = "AAAGGCTGGGGACCCTGATCTAAATACACCAATAAAAAGAAAAAGATTGTAAGATTGGATTTTAAAAGACCTGACTCTATACTGACCACAAAAAAACC"
    
    R_al,Q_al = smith_waterman(R,Q)
    t_end = time.time()
    print("Total execution time in (ms): ", 1000 * (t_end - t_start))

    print(R_al)
    print(Q_al)