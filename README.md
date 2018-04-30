# motifSearch

Author:  Dario Ghersi (dghersi[at]unomaha.edu)

A GUI program written in ```Python/QT4``` that calculates gene family distributions and clonotypes for CDR regions that match a given motif.
The program supports basic regular expression pattern matching and requires a tab-separated input file with deep sequencing data.

## Input data

See "sampleAlpha.tsv" for an example of how the data should look like:

1. Sample ID (e.g., "E1603")
2. Epitope (e.g., "BM")
3. Infection time point (e.g., "acute")
4. Type (e.g., "persistent")
5. CDR3 DNA sequence
6. CDR3 Amino acid sequence
7. Number of reads
8. CDR3 length (calculated as the length of the string minus two)

## Running the program

An executable version of the program for Mac OS X is available as "motifSearchMac". To run it, right-click on the file and select "Open With" -> "Terminal". In alternative, the program can be run by issuing the following command:

```
python motifSearch.py
```

The program requires the ```PyQT4``` package.

Load the input data by clicking on the "Load data", fill out all the entries and press the "Go!" button.

## Results

The "Matching Sequences" panel contains the DNA sequence, the amino acid sequence, the gene family/genes of the sequences matching the amino acid motif.

The "Gene Families and Usage" panel contains summative statistics on these sequences, including the proportion of gene family/genes with the motif (calculated over the clonotypes that match all search criteria but regardless of the presence/absence of the motif), by clonotype and by frequency.
