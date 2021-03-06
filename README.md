# Smith_Waterman_Alignment
Smith-Waterman Alignment of two input sequences 

This R-script re-creates the Smith-Waterman Alignment as described here (https://en.wikipedia.org/wiki/Smith–Waterman_algorithm). 

The script may be run as follows: 
`Rscript --vanilla hw1.R input.txt blosum62.txt`

Necessary inputs include: 

1. inputFile (example: input.txt): input sequence with 2 lines corresponding to 2 sequences to be aligned
2. scoreFile (example: blosum62.txt): score file with a matrix of each letter and the baseline score associated with each combination of letters.

Two additional optional inputs may be included for the opening gap and gap extension terms. Defaults for these are: openGap=-2, extGap=-1.

Outputs will include: 
1. The best alignment score. (Printed as ouput and also included in the filenames of the two outputs below:)
2. The full score matrix. (Sample output corresponding to the example input file is `AlignmentScores_283.output_scoreFile.tsv`)
3. The full alignment. (Sample output corresponding to the example input file is `AlignmentScore_283.output_Alignment.tsv`)

This page may be cloned using: 

```library(devtools)```
```install_github("mzekavat/Smith_Waterman_Alignment")```




