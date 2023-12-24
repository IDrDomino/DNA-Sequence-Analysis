# DNA-Sequence-Analysis-

## Overview

This Python program is designed for the analysis of DNA sequences. It leverages the capabilities of the Python programming language along with popular data science libraries such as NumPy and Pandas for efficient data manipulation and analysis.

## Environment
The code is intended to run in a Python 3 environment and is configured to work with the Kaggle/python Docker image. This environment comes equipped with various analytics libraries that enhance the functionality of the program. Notable dependencies include:

- **NumPy:** Used for efficient handling of arrays and mathematical operations.
- **Pandas:** Facilitates data processing and CSV file input/output.

```python
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import os
for dirname, _, filenames in os.walk('/kaggle/input'):
    for filename in filenames:
        print(os.path.join(dirname, filename))
```

## Bioinformatics Libraries

The following bioinformatics libraries are utilized for various tasks:

- BioPython: A comprehensive library for biological computation, providing tools for the manipulation and analysis of biological data. It includes modules for handling sequences, alignments, BLAST searches, and more.

- nglview: A library for molecular visualization that enhances the visual representation of molecular structures.

- matplotlib: A widely-used plotting library for creating static, animated, and interactive visualizations in Python.

- colorama: A library for adding colored output to terminal text, enhancing the readability of console messages.
  
- seaborn: A statistical data visualization library based on matplotlib, providing a high-level interface for drawing attractive and informative statistical graphics.

- pandas: A powerful data manipulation and analysis library that facilitates working with structured data.

```python
!pip install nglview
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.pairwise2 import format_alignment 
from Bio.Align import MultipleSeqAlignment
from Bio.Blast import NCBIWWW
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import CodonTable
import nglview as nv
import matplotlib.pyplot as plt
%matplotlib inline
from colorama import Back, Style, Fore
import seaborn as sns
import pandas as pd
```

## What is DNA?

Deoxyribonucleic acid, abbreviated as DNA, serves as the genetic code dictating the traits of living organisms. The structural components of DNA consist of nucleotides, each composed of a sugar and a phosphate molecule forming the DNA's 'backbone.' Additionally, DNA contains four organic bases: adenine (A), guanine (G), cytosine (C), and thymine (T).

## Sequence Representation Utility
This Python script includes functions for representing DNA and protein sequences with color-coded visualizations. The script provides two main functions:

```python
def ten_nucleotide_seq(genome):
    genes = []
    for ix, char in enumerate(genome):
        if ix != 0 and ix%10 == 0:
            genes.append(' ')
        genes.append(char)
    return ''.join(genes)

# color code to represent genome sequences
nu_clr_switcher = {
    # standard color-codes
    'A': Back.GREEN,
    'C': Back.YELLOW,
    'G': Back.RED,
    'T': Back.BLUE,
    ' ': Style.RESET_ALL
}
protein_clr_switcher = {
    # color-code by proteinfamily's polarity
    'A': Back.BLUE,
    'V': Back.BLUE,
    'I': Back.BLUE,
    'L': Back.BLUE,
    'M': Back.BLUE,
    'F': Back.BLUE,
    'Y': Back.CYAN,
    'W': Back.BLUE,
    'H': Back.CYAN,
    'R': Back.RED,
    'K': Back.RED,
    'N': Back.GREEN,
    'Q': Back.GREEN,
    'E': Back.MAGENTA,
    'D': Back.MAGENTA,
    'S': Back.GREEN,
    'T': Back.GREEN,
    'G': Back.YELLOW,
    'P': Back.YELLOW,
    'C': Back.BLUE,
    ' ': Style.RESET_ALL
}
def seq_repr(genome_str, strand ='dna'):
    if strand == 'dna':
        genome_str = ten_nucleotide_seq(genome=genome_str)
        line_break_cntr = 0
        for i in range(len(genome_str)):
            if genome_str[i] == ' ':
                line_break_cntr += 1
                if line_break_cntr>0 and line_break_cntr%6==0:
                    text = "\n"
                else:
                    text = nu_clr_switcher[genome_str[i]] + genome_str[i]
            else:
                text = nu_clr_switcher[genome_str[i]] + genome_str[i]
            print(text, end="")
        Style.RESET_ALL
    if strand == 'protein':
        for i in range(len(genome_str)):
            if genome_str[i] in protein_clr_switcher:
                if genome_str[i] == 'S' and genome_str[i+1:i+4] == 'TOP':
                    text = Style.RESET_ALL + 'S'
                elif genome_str[i] == 'T' and genome_str[i-1] == 'S' and genome_str[i+1:i+3] == 'OP':
                    text = Style.RESET_ALL + 'T'
                elif genome_str[i] == 'P' and genome_str[i-3:i] == 'STO':
                    text = Style.RESET_ALL + 'P'
                else:
                    text = protein_clr_switcher[genome_str[i]] + genome_str[i]
            else:
                Style.RESET_ALL
                text = genome_str[i]
            print(text, end="")
```

## Visualizing the Nucleotides of the Covid19 DNA 👇

```python
print("COVID-19 genome: ")
seq_repr(ncov_dna[0:300])
```
## Color Codes
DNA Sequence (strand='dna'):

- A (Adenine): Green
- C (Cytosine): Yellow
- G (Guanine): Red
- T (Thymine): Blue
- Space: Reset
- Protein Sequence (strand='protein'):

Amino acids are color-coded based on their properties, providing visual distinctions between different types.

#### COVID-19 genome 👇:  

![image](https://github.com/IDrDomino/DNA-Sequence-Analysis-/assets/154571800/347fe3ba-139c-47f8-8c91-84fb3c4a805f)

Total Number of Nucleotides

```python
len(ncov_dna)
#Output: 29845
```

## Finding the Composition of each Nucleotide
This Python function, nucleotides_composition, is designed to analyze the nucleotide composition of a given DNA sequence. It calculates the percentage of occurrence for each nucleotide (A, C, G, T) in the sequence.

```python
def nucleotides_composition(seq):
    nucleotides = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for n in nucleotides:
        nucleotides[n] = seq.count(n)/len(seq)*100
    return nucleotides
```

```python
ndict=nucleotides_composition(ncov_dna)
ndict
```

Insights into the distribution of nucleotides within the COVID-19 DNA sequence

```
{'A': 29.81738984754565,
 'C': 18.36153459540962,
 'G': 19.604623890098843,
 'T': 32.125984251968504}
```

```python
import pandas as pd
ndf = pd.DataFrame.from_dict(ndict, orient ='index')
ndf = ndf.reset_index()
ndf = ndf.rename(columns={"index": "Nucleotide", 0: "Composition"})
```

