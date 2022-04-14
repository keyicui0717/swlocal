
### Usage: python3 swlocal.py -i <input file> -s <score file>
### Example: python3 hw1.py -i sample-input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import numpy as np
import argparse


# Read in command-line arguments
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--similarity', help='similarity matrix file', required=True)
parser.add_argument('-og', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-eg', '--extgap', help='extension gap', required=False, default=-1)
parser.add_argument('-o', '--output', help='output gap', required=False, default="sw-alignment-output.txt")
args = parser.parse_args()


OPEN_PENALTY = int(args.opengap)
EXTENSION_PENALTY = int(args.extgap)


def get_sim_mat_dict(similiarity_mat_file):

    with open(similiarity_mat_file) as f:
        mat_header = f.readline().split()
        sim_mat_dict = {}
        for row_index, row in enumerate(f):
            row = row.strip()
            if not row:
                continue
            row = row.split()[1:]
            
            r = {}
            for column_index, score in enumerate(row):
                r[mat_header[column_index]] = score
            sim_mat_dict[mat_header[row_index]] = r
    return sim_mat_dict

def get_mat(seqx, seqy, sim_mat_dict):
    column_count = len(seqx) + 1
    row_count = len(seqy) + 1
    sw_mat = np.zeros((row_count, column_count), int)
    tb_mat = np.zeros((row_count, column_count, 2), int)


    sw_mat_max = [0, (None, None)] 
    for i in range(1,row_count):
        for j in range(1, column_count):
            elem1 = seqy[i-1]
            elem2 = seqx[j-1]
            similarity_score = int(sim_mat_dict[elem1][elem2])
            sw_mat[i, j], tb_mat[i, j] = get_sw_score(similarity_score, sw_mat, i, j)
            sw_mat_max = [sw_mat[i, j], [i, j]] if sw_mat_max[0] < sw_mat[i,j] else sw_mat_max

    return sw_mat, tb_mat, sw_mat_max

def get_sw_score(similarity_score, sw_mat, i, j):
    up = [(sw_mat[i-len][j] + gap_penalty(len)) for len in range(1, i+1)]
    left = [(sw_mat[i][j-len] + gap_penalty(len)) for len in range(1, j+1)]

    from_diagonal = similarity_score + sw_mat[i-1][j-1]

    max_from_up = max(up)
    for count, score in enumerate(up):
        if score == max_from_up:
            gap_from_up = count + 1 

    max_from_left = max(left)
    for count, score in enumerate(left):
        if score == max_from_left:
            gap_from_left = count + 1 

    scores = (from_diagonal, max_from_up, max_from_left)

    if max(scores) < 0:
        sw_score = 0
        return sw_score, [0, 0] #list is [direction, distance]
    else:
        sw_score = max(scores)
    
    tb_direction = np.argmax(scores)+1
    if tb_direction == 1:
        return sw_score, [1, 1]
    elif tb_direction == 2:
        return sw_score, [2, gap_from_up]
    elif tb_direction == 3:
        return sw_score, [3, gap_from_left]
    

def gap_penalty(length):
    """ Returns the cost of a gap dependent on the length of the gap
        and cost for opening a gap versus extending a gap  """
    if length == 1:
        return OPEN_PENALTY
    elif length > 1:
        return OPEN_PENALTY + EXTENSION_PENALTY * (length-1)

def get_alignment(initial_i, initial_j, sw_mat, tb_mat, seqx, seqy):
    
    # Traceback 
    i = initial_i
    j = initial_j
    
    tb_seqx = ""
    tb_seqy = ""

    while sw_mat[i, j] != 0:
        tb_direction = tb_mat[i,j][0]
        tb_gap = tb_mat[i,j][1]

        if tb_direction == 1:
            tb_seqx += (seqx[j-1])
            tb_seqy += (seqy[i-1])
            i, j = i-1, j-1

        #if direction was from up, it means gap in seqx
        elif tb_direction == 2:
            tb_seqx += ("-" * tb_gap)
            tb_seqy += "".join(reversed(seqy[i-tb_gap:i]))
            i,j = i-tb_gap, j
        
        #if direction was form left, it means a gap in seqy
        elif tb_direction == 3:
            tb_seqx += "".join(reversed(seqx[j-tb_gap:j]))
            tb_seqy += ("-" * tb_gap)
            i, j = i, j-tb_gap

    
    # The subsequences in each sequence that are aligned obtained via traceback
    tb_seqx = ''.join(reversed(tb_seqx))
    tb_seqy = ''.join(reversed(tb_seqy))
    

    # join the rest of the sequence surrounding locally aligned subsequences 

    align_seqx = seqx[ :j] + "(" + tb_seqx + ")" + seqx[initial_j:]
    align_seqy = seqy[ :i] + "(" + tb_seqy + ")" + seqy[initial_i: ]

    if j > i : # if seq1 aligns at a further position relative to the position in seq2 where alignment starts
        pad = (j-i) * " "
        align_seqy = pad + align_seqy  #pad beginning of second sequence with spaces
 
    elif i > j: # if seq2 aligns at a further position relative to the position in seq1 where alignment starts
        pad = (i-j) * " "
        align_seqx = pad + align_seqx #pad beginning of first sequence with spaces

    # Indidate which residues match
    match = ""
    for n in range(min(len(align_seqx), len(align_seqy))):
        if align_seqx[n]!= align_seqy[n] or align_seqx[n] in ["(", ")"]:
            match += " "
        elif align_seqx[n] == align_seqy[n]:
            match += "|"

    return align_seqx, align_seqy, match


def runSW(input_file, similiarity_mat_file):
    with open(input_file) as f:
        seqx = next(f).strip() #sequence 1
        seqy = next(f).strip() #sequence 2

    sim_mat_dict = get_sim_mat_dict(similiarity_mat_file)
    sw_mat, tb_mat, sw_mat_max = get_mat(seqx, seqy, sim_mat_dict)
    
    alignment_score = sw_mat_max[0]
    tb_inital_i = sw_mat_max[1][0]
    tb_initial_j = sw_mat_max[1][1]

    align_seqx, align_seqy, match = get_alignment(tb_inital_i, tb_initial_j, sw_mat, tb_mat, seqx, seqy)

    print(f"alignment score is: {alignment_score}")
    print(f"Alignment Results:\n{align_seqx}\n{match}\n{align_seqy}")

    with open(args.output,"a") as output:
        output.write("-----------\n|Sequences|\n-----------\n")
        output.write(f"sequence1\n{seqx}\nsequence2\n{seqy}\n")
        output.write("--------------\n|Score Matrix|\n--------------\n")
        output.write("\t\t" + "\t".join([elem for elem in seqx]) + "\n")
        for count, row in enumerate(sw_mat):
            if count == 0:
                output.write("\t" + "\t".join([str(score)for score in row]) + "\n")
            else:
                output.write(seqy[count-1] + "\t" + "\t".join([str(score)for score in row]) + "\n")
        output.write("----------------------\n|Best Local Alignment|\n----------------------\n")
        output.write(f"Alignment Score:{alignment_score}\n")
        output.write(f"Alignment Results:\n{align_seqx}\n{match}\n{align_seqy}\n")
        
if __name__ == "__main__":
    runSW(args.input, args.similarity)