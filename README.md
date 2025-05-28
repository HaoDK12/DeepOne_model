# DeepOne

A deep-learning model for SpCas9-induced +1bp frequency in cells

## Introduction

DeepOne is an online computational tool designed to predict templated 1-bp insertion frequencies in CRISPR/Cas9-mediated genome editing. Previous studies have domonstrated that CRISPR editing precision correlates 
with not only overall editing activity but also a strong preference for single-nucleotide homologous insertion. These insertions predominantly occur at the 17th nucleotide position (N17), four bases upstream of the PAM 
sequence, where the inserted nucleotide is typically a duplicate of N17 in the target sequence. Hopefully, this precise and templated repair is considered to be desirable strategy, particularly for correcting pathogenic 
alleles and performing gene knockouts. Here, we leveraged templated 1-bp insertion frequency data from over 15,000 gRNAs in HEK293 cells, and developed a deep learning-based tool, specifically optimized for this task. 
Considering that insertion frequency varies across different cellular contexts, we further fine-tuned the model using data from Allen et al. for six additional cell lines: K562 human chronic myelogenous leukemia cells, 
human induced pluripotent stem cells (iPSCs), mouse embryonic stem cells (mESCs), Chinese hamster ovary (CHO) cells, human retinal pigment epithelial cells (RPE-1), and leukemic near-haploid (HAP1) cells. In independent 
benchmarking, DeepOne outperformed existing model, domonstrating state-of-the-art accuracy in predicting templated 1-bp insertion probabilities acorss multiple cell types. A user-friendly computational portal is available
at https://dreamdb.biomed.au.dk/DeepOne/home.

## Requirements
The scripts are written in Python 3.7.16 and run on LINUX. The versions of python packages which we used are, specifically:
``` 
numpy_1.21.5 pandas_1.3.5  tensorflow_2.11.0
keras_2.11.0   h5py_3.7.0  viennarna_2.3.3
```

## Installation Guide
Clone this github repository, then set up your environment to import the DeepOne-model.py script in however is most convenient for you. In python, for instance, you may use the following at the top of your script to import DeepOne.
```
import sys
sys.path.append('/directory/to/local/DeepOne-model/repo/clone/')
import DeepOne-model
```

## Run DeepOne with command line tool
Run the DeepOne-model command line interface (CLI) with the command in the terminal:
```
python DeepOne-model.py --input_seq TTATCTTCGCTATCACCTCCGCCGGGGTCACCCATTAT --cell_line HEK293 --out_path results.tsv --prefix sample_
```
Note: Input seq is a genomic sequence limited from 31 to 2,000 bp, ensuring all spaces, line breaks, and numbers are removed. Supported cell types are ['mESC', 'CHO', 'HEK293', 'IPSC', 'K562', 'HAP1', 'RPE-1']. Given that DeepOne-HEK model offer superior predictive performance acorss multiple cell type in our analysis, we recommend using HEK293 cell if your cell type of interest is not listed here. Alternatively, several similar cells, such as Human embryonic stem cells and mESCs, are considerable.

##
