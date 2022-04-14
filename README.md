# Smith-Waterman Local Alignment

## The Smith-Waterman algorithm preforms local alignment between two sequences. It is a variaiton of a dynamic programming algorithm and finds the optimal local alignment with respect to the specified scoring system.



## Installation
```bash
pip install git+https://github.com/keyicui0717/swlocal.git
```
## Requirements
* numpy
  
## Usage
Run program from command-line.
```bash
python3 swlocal.py -i <input_file> -s <similarity_matrix_file>
```
## Example
```bash
python3 swlocal.py -i sample-input1.txt -s blosum62.txt
```
Example input
```txt
DEMKASEDLKKHGATVLTAL
VSAHGATVLNKLG
```
Example output
```txt
-----------
|Sequences|
-----------
sequence1
DEMKASEDLKKHGATVLTAL
sequence2
VSAHGATVLNKLG
--------------
|Score Matrix|
--------------
		D	E	M	K	A	S	E	D	L	K	K	H	G	A	T	V	L	T	A	L
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
V	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	4	2	1	0	1
S	0	0	0	0	1	1	4	2	1	0	1	0	0	0	1	1	2	2	3	2	0
A	0	0	0	0	0	5	3	3	1	0	0	0	0	0	4	2	1	1	2	7	5
H	0	0	0	0	0	3	4	3	2	0	0	0	8	6	5	4	3	2	1	5	4
G	0	0	0	0	0	2	3	2	2	0	0	0	6	14	12	11	10	9	8	7	6
A	0	0	0	0	0	4	3	2	0	1	0	0	5	12	18	16	15	14	13	12	11
T	0	0	0	0	0	2	5	3	2	1	0	0	4	11	16	23	21	20	19	18	17
V	0	0	0	1	0	1	3	3	1	3	1	0	3	10	15	21	27	25	24	23	22
L	0	0	0	2	0	0	2	1	0	5	3	2	2	9	14	20	25	31	29	28	27
N	0	1	0	0	2	0	1	2	2	3	5	3	3	8	13	19	24	29	31	29	28
K	0	0	2	0	5	3	2	2	1	2	8	10	8	7	12	18	23	28	29	30	28
L	0	0	0	4	3	4	2	1	0	5	6	8	7	6	11	17	22	27	28	28	34
G	0	0	0	2	2	3	4	2	1	3	5	7	6	13	11	16	21	26	27	28	32
----------------------
|Best Local Alignment|
----------------------
Alignment Score:34
Alignment Results:
DEMKASEDLKK(HGATVLTAL)
            ||||||  | 
        VSA(HGATVLNKL)G
```