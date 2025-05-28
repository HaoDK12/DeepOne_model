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
