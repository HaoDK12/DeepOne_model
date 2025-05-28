# DeepOne

A deep-learning model for SpCas9-induced +1bp frequency in cells

## Introduction

DeepOne is an online and offline computational tool designed to predict templated 1-bp insertion frequencies in CRISPR/Cas9-mediated genome editing. Based on previous findings<sup>1,2,3</sup>, such insertions—especially at the 17th position upstream of the PAM—are frequent and biologically meaningful. Leveraging data from over 15,000 gRNAs, we trained a deep learning model in HEK293 cells, and fine-tuned it for six other cell lines<sup>4</sup>.

> DeepOne achieves state-of-the-art accuracy across various cellular contexts and is freely available at [dreamdb.biomed.au.dk/DeepOne](https://dreamdb.biomed.au.dk/DeepOne/home).

## Requirements
Test with
``` 
numpy_1.21.5 pandas_1.3.5  tensorflow_2.11.0
keras_2.11.0   h5py_3.7.0  viennarna_2.3.3
```

## Installation
```bash
git clone git clone https://github.com/HaoDK12/DeepOne_model.git
cd DeepOne_model
```
Set environment as needed. You can import and use the model in scripts or run it via CLI.

## Run DeepOne via command line tool
Run the DeepOne-model command line interface (CLI) with the command in the terminal:
```
python DeepOne-model.py --input_seq TTATCTTCGCTATCACCTCCGCCGGGGTCACCCATTAT --cell_line HEK293 --out_path results.tsv --prefix sample_
```
| Argument       | Type  | Required | Description                                                                                          |
| -------------- | ----- | -------- | ---------------------------------------------------------------------------------------------------- |
| `--input_seq`  | `str` | Yes    | Genomic DNA sequence (31–2000 bp) without spaces, line breaks, or numbers.                           |
| `--cell_line`  | `str` | Yes    | Cell line used for prediction. Supported: `HEK293`, `CHO`, `HAP1`, `iPSC`, `K562`, `mESCs`, `RPE-1`. |
| `--out_path`   | `str` | Yes    | Output file path for saving the prediction results (e.g., `results.tsv`).                            |
| `--prefix`     | `str` | No     | Optional prefix for the guide ID column (default: none).                                             |
| `--help`, `-h` | flag  | No     | Show help message and exit.                                                                          |

Note: Supported cell types are ['mESC', 'CHO', 'HEK293', 'IPSC', 'K562', 'HAP1', 'RPE-1']. Given that DeepOne-HEK model offer superior predictive performance acorss multiple cell type in our analysis, we recommend using HEK293 cell if your cell type of interest is not listed here. Alternatively, several similar cells, such as Human embryonic stem cells and mESCs, are considerable.

## Contact
We greatly appreciate your feedback. If bug reports or suggestions, Please contact us (au735018@uni.au.dk).

## Cite
<sup>1</sup> Chakrabarti AM, Henser-Brownhill T, Monserrat J, Poetsch AR, Luscombe NM, Scaffidi P. Target-Specific Precision of CRISPR-Mediated Genome Editing. Mol Cell. 2019 Feb 21;73(4):699-713.e6.  
<sup>2</sup> Shen MW, Arbab M, Hsu JY, Worstell D, Culbertson SJ, Krabbe O, Cassa CA, Liu DR, Gifford DK, Sherwood RI. Predictable and precise template-free CRISPR editing of pathogenic variants. Nature. 2018 Nov; 563(7733):646-651.  
<sup>3</sup> Shou J, Li J, Liu Y, Wu Q. Precise and Predictable CRISPR Chromosomal Rearrangements Reveal Principles of Cas9-Mediated Nucleotide Insertion. Mol Cell. 2018 Aug 16; 71(4):498-509.e4.  
<sup>4</sup> Chen W, McKenna A, Schreiber J, Haeussler M, Yin Y, Agarwal V, Noble WS, Shendure J. Massively parallel profiling and predictive modeling of the outcomes of CRISPR/Cas9-mediated double-strand break repair. Nucleic Acids Res. 2019 Sep 5; 47(15):7989-8003.

If you are using DeepOne in your publication, please cite:  
Hao Yuan, Xiaoguang Pan, Menachem Viktor Khamo Sarusie, Janos Hasko, Huixin Xu, Chunping Song, Julie Lund Petersen, Trine Skov Petersen, Soren Tvorup Christensen, Lars Allan Larsen, Lin Lin, Yonglun Luo#, Precise and efficient CRISPR gene 
editing by 1bp insertion with deep-learning gRNAs. 2025 (Manuscript under revision)
