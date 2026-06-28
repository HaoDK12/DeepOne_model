# DeepOne

A deep-learning model for SpCas9-induced +1bp frequency in cells

## Introduction

DeepOne is an online and offline computational tool designed to predict templated 1-bp insertion frequencies in CRISPR/Cas9-mediated genome editing. Based on previous findings<sup>1,2,3</sup>, such insertions—especially at the 17th position upstream of the PAM—are frequent and biologically meaningful. Leveraging data from over 15,000 gRNAs, we trained a deep learning model in HEK293 cells, and fine-tuned it for six other cell lines<sup>4</sup>. To broaden the applicability of DeepOne, we extended it to SpCas9-NG and SpG 1-bp prediction tasks.

> DeepOne achieves state-of-the-art accuracy across various cellular contexts and SpCas9 enzymes and is freely available at [dreamdb.biomed.au.dk/DeepOne](https://dreamdb.biomed.au.dk/DeepOne/home).

## Requirements
Test with
``` 
numpy_1.21.5 pandas_1.3.5  tensorflow_2.11.0
keras_2.11.0   h5py_3.7.0  viennarna_2.3.3
```

## Installation
```bash
git clone https://github.com/HaoDK12/DeepOne_model.git
cd DeepOne_model
```
Set environment as needed. You can import and use the model in scripts or run it via command line interface (CLI).

## Run DeepOne via command line tool
Run the DeepOne-model CLI with the command in the terminal:
```
python DeepOne-model.py --input_seq TTATCTTCGCTATCACCTCCGCCGGGGTCACCCATTAT --cell_line HEK293 --variant --out_path results.tsv --prefix sample_
```
| Argument       | Type  | Required | Description                                                                                          |
| -------------- | ----- | -------- | ---------------------------------------------------------------------------------------------------- |
| `--input_seq`  | `str` | Yes    | Genomic DNA sequence (31–2000 bp) without spaces, line breaks, or numbers.                           |
| `--cell_line`  | `str` | Yes    | Cell line used for prediction. Supported: `HEK293`, `CHO`, `HAP1`, `iPSC`, `K562`, `mESCs`, `RPE-1`. |
| `--variant`  | `str` | No   | Optional model prediction for cas9 variants. Supported: `None`,`SpCas9-NG`,`SpG`, `Both`. |
| `--out`   | `str` | Yes    | Output file path for saving the prediction results (e.g., `Output.tsv`).                            |
| `--prefix`     | `str` | No     | Optional prefix for the guide ID column (default: none).                                             |
| `--help`, `-h` | flag  | No     | Show help message and exit.                                                                          |

Note: 
1) Supported cell types are ``['mESC', 'CHO', 'HEK293', 'IPSC', 'K562', 'HAP1', 'RPE-1']``. Given that DeepOne-HEK model offers superior predictive performance across multiple cell types in our analysis, we recommend using HEK293 cell if your cell type of interest is not listed here. Alternatively, several similar cells, such as Human embryonic stem cells and mESCs, may also be considered. 
2) Supported cas9 variants are provided as optional selections, and the corresponding variant prediction will be included in the output if chosen.
3) Prediction outcome might vary slightly because of differences in the version or installment of the Viennarna package, which can affect the calculated input parameter delta-Gb.

## Advanced Workflow: Joint prediction pipeline of 1-bp insertion frequencies (DeepOne + inDelphi)

To select the overall optimal sgRNA candidates for frame-restoration or precise knockouts, we provide an automated joint prediction pipeline. This workflow screens both strands of a long genomic region for valid NGG PAM sites, computes absolute $+1$ bp duplicated insertion frequencies using **DeepOne**, profiles indel spectrum landscapes via **inDelphi**, and ultimately ranks candidates based on the custom +1bp ins precision index:
```math
$$\text{Ins1bp\_vs\_Frameshift\_Ratio} = \frac{\text{inDelphi 1-bp Insertion Frequency}}{\text{inDelphi Total Frameshift Frequency}}$$
```
This helps filter out targets where $+1$ bp insertions are diluted by undesirable indels.

### Decoupled Environment Setup
Because `DeepOne` depends on frameworks (Python 3 / TensorFlow 2.11) and `inDelphi` runs on legacy dependencies (Python 2.7), the pipeline isolates these tasks seamlessly via two sub-process bridging.

1. Clone both repositories into your local system:
   - DeepOne: `https://github.com/HaoDK12/DeepOne_model.git`
   - inDelphi: `https://github.com/maxwshen/inDelphi-model.git`
2. Ensure you have corresponding Python 3 and Python 2.7 environments activated or mapped.

### Configuration (`config.ini`)
Parameters are managed externally via `config.ini`. Create this file in the same workspace directory as `run_pipeline.py` and specify your paths and targets:

```ini
[SEQUENCE_AND_CELL]
SAMPLE_GENOMIC_SEQUENCE = GCTATCACCTCCGCCGGGGTCACCCATTATCTGGCGGCCCCCCGGAAAGGGGGGGGGGGGGG
CELL_LINE = HEK293

[DEEPONE_ENV]
PYTHON3_EXEC = python
DEEPONE_DIR = /path/to/local/DeepOne_model/

[INDELPHI_ENV]
PYTHON2_EXEC = /path/to/py27/bin/python
INDELPHI_DIR = /path/to/local/inDelphi-model-master/
```
| Argument       | Required | Description                                                                                          |
| -------------- | -------- | ---------------------------------------------------------------------------------------------------- |
| `--SAMPLE_GENOMIC_SEQUENCE`  | str (>= 60 bp genomic continuous DNA)    | Target window containing flanking sequences around PAM sites. |
| `--CELL_LINE`  | "HEK293, mESC, K562"    | Limited to the intersection supported by both pre-trained models. |
| `--PYTHON3_EXEC / PYTHON2_EXEC`  | environment paths  | Paths to the respective Python 3 and Python 2.7 runtime binaries. |
| `--DEEPONE_DIR / INDELPHI_DIR`  | absolute local directories  | Root path of cloned Git repositories. |


## Running the Joint Pipeline
Execute the master script in your terminal using Python 3:
```
python run_pipeline_master.py
```
## Output
The workflow evaluates targets on both strands and generates a tab-separated table file `./Final_Joint_Predictions.tsv` containing fused rows (ID, Target, Strand, PAM, GC%, DeepOne Score, inDelphi 1bp Ins%, inDelphi Frameshift%). The entries are hierarchically sorted in descending order based on the `Ins1bp_vs_Frameshift_Ratio` to instantly highlight elite target sites.


## Contact
We greatly appreciate your feedback. If bug reports or suggestions, Please contact us (yuanhao971@gmail.com).

## Cite
<sup>1</sup> Chakrabarti AM, Henser-Brownhill T, Monserrat J, Poetsch AR, Luscombe NM, Scaffidi P. Target-Specific Precision of CRISPR-Mediated Genome Editing. Mol Cell. 2019 Feb 21;73(4):699-713.e6.  
<sup>2</sup> Shen MW, Arbab M, Hsu JY, Worstell D, Culbertson SJ, Krabbe O, Cassa CA, Liu DR, Gifford DK, Sherwood RI. Predictable and precise template-free CRISPR editing of pathogenic variants. Nature. 2018 Nov; 563(7733):646-651.  
<sup>3</sup> Shou J, Li J, Liu Y, Wu Q. Precise and Predictable CRISPR Chromosomal Rearrangements Reveal Principles of Cas9-Mediated Nucleotide Insertion. Mol Cell. 2018 Aug 16; 71(4):498-509.e4.  
<sup>4</sup> Chen W, McKenna A, Schreiber J, Haeussler M, Yin Y, Agarwal V, Noble WS, Shendure J. Massively parallel profiling and predictive modeling of the outcomes of CRISPR/Cas9-mediated double-strand break repair. Nucleic Acids Res. 2019 Sep 5; 47(15):7989-8003.

If you are using DeepOne in your publication, please cite:  
Hao Yuan, Xiaoguang Pan, Menachem Viktor Khamo Sarusie, Janos Hasko, Huixin Xu, Chunping Song, Julie Lund Petersen, Trine Skov Petersen, Soren Tvorup Christensen, Lars Allan Larsen, Lin Lin, Yonglun Luo#, Precise and efficient CRISPR gene 
editing by 1bp insertion with deep-learning gRNAs. 2025 (Manuscript under revision)
