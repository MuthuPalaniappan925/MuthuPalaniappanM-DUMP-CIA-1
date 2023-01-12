def cia_affine_model(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend):
    
    m, n = len(seq1), len(seq2)
    
    #1-> Step 1 To create a scoring matrix
    # Create the scoring matrix
    
    score_matrix = [[0] * (n + 1) for i in range(m + 1)]
    
    #2-> Step 2 To create a tracback matrix
    # Create the traceback matrix
    
    traceback_matrix = [[0] * (n + 1) for i in range(m + 1)]
    
    #assigning the maximum score as zero and position array.

    max_score = 0
    max_pos = None
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate the match/mismatch score
            if seq1[i - 1] == seq2[j - 1]:
                diag_score = score_matrix[i - 1][j - 1] + match_score
            else:
                diag_score = score_matrix[i - 1][j - 1] + mismatch_score
            # Calculate the gap scores
            gap1 = score_matrix[i - 1][j] + gap_open + gap_extend
            gap2 = score_matrix[i][j - 1] + gap_open + gap_extend
            gap3 = score_matrix[i - 1][j - 1] + gap_open + gap_extend
            # Find the maximum score
            score = max(0, diag_score, gap1, gap2, gap3)
            score_matrix[i][j] = score
            if score == 0:
                traceback_matrix[i][j] = 0
            elif score == gap1:
                traceback_matrix[i][j] = 1
            elif score == gap2:
                traceback_matrix[i][j] = 2
            elif score == diag_score:
                traceback_matrix[i][j] = 3
            elif score == gap3:
                traceback_matrix[i][j] = 4
            if score > max_score:
                max_score = score
                max_pos = (i, j)

    # Traceback to find the alignment
    align1, align2 = "", ""
    i, j = max_pos
    while traceback_matrix[i][j] != 0:
        if traceback_matrix[i][j] == 3:
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 1:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
        elif traceback_matrix[i][j] == 2:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
        elif traceback_matrix[i][j] == 4:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
    #printing out the traceback matrix
    print(f"Sequence 1 : {x}\nSequence 2 : {y}")
    print("\n")
    print("Traceback Matrix\n")
    print(traceback_matrix)
    print("\n")
    print("Score Matrix\n")
    print(score_matrix)
    print("\n")
    print("Alligned Sequence\n")
    print(align1)
    print(align2)
    print("\n")
    print(f"Match Score : {match_score}\nMis-Match Score : {mismatch_score}\nGap Open : {gap_open}\nGap Extend : {gap_extend}")
    return align1,align2

x = "AAAACCCCGGGGTTTT"
y = "TAAAACCCCGGGGTTT"

a1,a2 = cia_affine_model(x, y, 5, -1, -5, -10)
