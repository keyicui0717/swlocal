# Smith-Waterman Local Alignment

## The Smith-Waterman algorithm preforms local alignment between two sequences. It is a variaiton of a dynamic programming algorithm finds the optial local alignment with respect to the specified scoring system.



## Installation
```bash
pip install git+https://github.com/keyicui0717/swlocal.git
```

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
MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG
VLKCWGPMEADYATHGGLVLTRLFTEHPETLKLFPKFAGIAHGDLAGDAGVSAHGATVLNKLGDLLKARGAHAALLKPLSSSHATKHKIPIINFKLIAEVIGKVMEEKAG
```
Example output
```txt




