# Direct Inference of Haplotypes from Sequencing Data

## 1. Introduction

This guide is intended to help readers utilize the haplotype inference tool ***DIHap***. Unlike previous algorithms, ***DIHap*** can simultaneously infer the sequencing error profile and haplotypes, thereby reducing accumulated errors. The accompanying paper demonstrates its robustness against higher sequencing errors and its capability to infer haplotypes in polyploid cases.

In bulk sequencing data, the reads lack cell type labels. The advantages of ***DIHap*** will assist researchers in studying the correlation between single nucleotide variations (SNVs), enabling multivariate association studies related to disease symptoms. Additionally, for certain animals and crops, such as fish and potatoes, ***DIHap*** can help identify relevant variations related to feed and cultivation. Thus, the algorithm has broad applicability in both research and practical applications.

***DIHap*** is run under the conda environment with python code.

## 2. Steps for installation

1. Follow the https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html to download the conda environment on your Linux or MacOS operation system. Although, the Windows system supports the conda, it does not work with some package related with biological package.

2. Down load the zipped package named "[DIHap_sandbox.tar.gz](https://drive.google.com/drive/folders/1jAqb9CF_JGgxz27k0paoil-e3OnzPami?usp=share_link)" from https://drive.google.com/drive/folders/1jAqb9CF_JGgxz27k0paoil-e3OnzPami?usp=share_link.

3. Upload and unzip the "samtools.rar" to the "envs" folder in your server or own computer. One example after Step 1 is "\home\anaconda3\envs\".

4. If your unzip the file and name "samtools". Then, after you step in the conda environment, you should use the code ```conda activate samtools``` to activate the environment. The reference material is from the "conda-pack" , see https://conda.github.io/conda-pack/.

## 3. File content

1. "Reads_SNV_inference.py" is the script for the inference with no long range variation case. And a folder named "code_for_SNV_inference" is the supporting functions in the script, which should be put in the same directory of the code.
2. "Reads_SNV_inference_accelerate_complex_version.py" is the script for the inference with long range variation case.
3. "Toy example" folder is the folder to contain the example files.

## 4. Usage

Usage of "Reads_SNV_inference.py" or "Reads_SNV_inference_accelerate_complex_version.py": ``` python Reads_SNV_inference.py main_url_bam main_url_ref_con main_url_save num_reads proposed_k num_cores``` or``` python Reads_SNV_inference_accelerate_complex_version.py main_url_bam main_url_ref_con main_url_save num_reads proposed_k num_cores```  , where ```main_url_bam``` is the path of bam file of the reads, ```main_url_ref``` is the path of the reference path, ```main_url_save``` is the saving files path, ```num_reads``` is the number of the reads and ```proposed_k``` is the number of the haplotypes, and ```num_cores``` is the number of cores for parallel computing.

Toy Example: ```python Reads_SNV_inference.py  \Toy example\0_reads_bam.bam \Toy example\0_reads_reference.fasta \Toy example 50 2 5``` to infer the 2 haplotypes of one person with 50 reads using 5 cores. 



