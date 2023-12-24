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

Now, we converted the nucleotide composition dictionary (ndict) into a Pandas DataFrame (ndf). The DataFrame includes two columns: "Nucleotide" and "Composition," making it a structured and easy-to-read format for further analysis or visualization.

```python
import pandas as pd
ndf = pd.DataFrame.from_dict(ndict, orient ='index')
ndf = ndf.reset_index()
ndf = ndf.rename(columns={"index": "Nucleotide", 0: "Composition"})
```

Now using the Seaborn library to create a bar plot based on the nucleotide composition data stored in the Pandas DataFrame (ndf).

```python
ax = sns.barplot(x="Nucleotide", y="Composition", data=ndf)
```

![__results___19_0](https://github.com/IDrDomino/DNA-Sequence-Analysis-/assets/154571800/4fbf6b13-378c-461b-8934-3a2d99625905)

# Calculating GC-content of the DNA
In polymerase chain reaction (PCR) studies, the GC-content of short oligonucleotides referred to as primers is frequently employed to estimate their annealing temperature to the template DNA. A heightened GC-content signifies a comparatively elevated melting temperature, thereby indicating increased stability.

GC function from BioPython's SeqUtils module to calculate the GC percentage of the DNA sequence (ncov_dna).
```python
from Bio.SeqUtils import GC
print(f"GC% :{GC(ncov_dna)}")
```
```
GC% :37.96615848550846
```

# Tri-nucleotide compositions
Within the realm of bioinformatics, k-mers are subsequences of a specified length denoted as {\displaystyle k}. This concept is integral to computational genomics and sequence analysis, particularly in the field of genomics where k-mers consist of nucleotides (A, T, G, and C). The application of k-mers is widespread, encompassing tasks such as DNA sequence assembly, enhancement of heterologous gene expression, species identification in metagenomic samples, and the development of attenuated vaccines.

Defined a set of trinucleotides (trimers) and created a function, trimer_composition, to calculate the count of each trinucleotide in a given DNA sequence (genome). This function will return a dictionary where keys are the trimers, and values are their respective counts in the provided DNA sequence.

```python
# tri-nucleotide compositions
trimers = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "ATA", "ATC", "ATG", "CAA", 
           "CAC", "CAG", "CCA","CCC","CCG","CGA","CGC","CTA","CTC","GAA","GAC","GCA","GCC","GGA","GTA","TAA","TCA"]

def trimer_composition(genome):
    trimer_dict = dict()
    for trimer in trimers:
        trimer_dict[trimer] = genome.count(trimer)
    return trimer_dict
```

```
composition = trimer_composition(ncov_dna)
total_composition = sum(composition.values())
norm_freq = [count/total_composition for count in composition.values()]
print(composition)
print(total_composition)
print(norm_freq)
```

```
{'AAA': 641, 'AAC': 615, 'AAG': 575, 'AAT': 760, 'ACA': 754, 'ACC': 371, 'ACG': 164, 'ACT': 674, 'AGA': 567, 'AGC': 300, 'AGG': 328, 'ATA': 446, 'ATC': 339, 'ATG': 722, 'CAA': 702, 'CAC': 427, 'CAG': 437, 'CCA': 353, 'CCC': 102, 'CCG': 73, 'CGA': 95, 'CGC': 94, 'CTA': 555, 'CTC': 270, 'GAA': 531, 'GAC': 338, 'GCA': 371, 'GCC': 188, 'GGA': 280, 'GTA': 468, 'TAA': 717, 'TCA': 550}
13807
[0.046425726080973416, 0.04454262330701818, 0.04164554211631781, 0.05504454262330702, 0.054609980444701965, 0.026870428043745925, 0.011878032881871515, 0.04881581806330122, 0.041066125878177734, 0.02172810893025277, 0.02375606576374303, 0.032302455276309115, 0.02455276309118563, 0.05229231549214167, 0.05084377489679148, 0.030926341710726443, 0.03165061200840154, 0.02556674150793076, 0.007387557036285942, 0.005287173173028174, 0.006880567827913377, 0.006808140798145868, 0.040197001520967626, 0.019555298037227494, 0.038458752806547404, 0.024480336061418122, 0.026870428043745925, 0.013616281596291736, 0.020279568334902586, 0.03389584993119432, 0.05193018034330412, 0.03983486637213008]`
```

Now, converting the trinucleotide composition dictionary (composition) into a Pandas DataFrame (tri)
```python
tri = pd.DataFrame.from_dict(composition, orient ='index')
tri = tri.reset_index()
tri = tri.rename(columns={"index": "trimer", 0: "count"})
```

Now, using Pandas and styling with a bar chart and a background gradient to visualize the trinucleotide counts in descending order.

```python
r1 = tri.sort_values(by='count', ascending=False)
r1.style.bar(subset=["count"],color='#').background_gradient(cmap='Reds')
```

![image](https://github.com/IDrDomino/DNA-Sequence-Analysis-/assets/154571800/db5ff49b-684c-4c00-b857-9eecea7274fa)

Finally, creates a horizontal bar plot where trinucleotides are represented on the y-axis, and their respective counts are represented on the x-axis. The fig_dims variable determines the size of the plot.

```python
fig_dims = (12, 8)
fig, ax = plt.subplots(figsize=fig_dims)
sns.barplot(x="count", y="trimer", ax=ax, data=tri)
```
![__results___27_1](https://github.com/IDrDomino/DNA-Sequence-Analysis-/assets/154571800/96178874-a6a5-4fdc-bca6-b2f85f9c312e)



