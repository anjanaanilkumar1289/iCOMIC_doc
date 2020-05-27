## 3. Getting started



This guide will walk you through the steps necessary to understand, install, and use iCOMIC for carrying out analysis on your data.

  

#### 3.1. iCOMIC overview

iCOMIC is an open-source, stand-alone toolkit for genomic data analysis, characterized by a python based Graphical User Interface. The tool enables researchers with minimal programming expertise to draw consequential insights from DNA-Seq and RNA-Seq data.

#### 3.2. Install iCOMIC

Installation is easy as we provide a `requirements.txt` file comprising all the software dependencies. Once you clone the iCOMIC github repository, you can install all the associated dependencies using the command below. Every additional software requirement will be managed by the conda environment.

```

$ pip install requirements.txt

```

  

#### 3.3. Usage

1) Snakemake

	iCOMIC is embedded in Snakemake, a Python based workflow manager. Different tools integrated in iCOMIC are connected using Snakemake. Individual ‘Rules’ corresponding to each tool form the building units, which describes how the desired output is obtained from the input. Rules consist of information about the input and output files and wrapper script or shell command. Tools without wrapper scripts are configured separately and shell command is used for their execution. According to the choice of tools made by the user, corresponding rules are combined in a ‘Snakefile’ to generate target output. All the input information and parameters corresponding to each tool is specified in a configuration file, ‘config file’. ‘Rules’ are predefined and are made available together with the iCOMIC package. All the other files are generated on the flow according to the user inputs and are updated accordingly. 

2) PyQt5 GUI

	iCOMIC is characterized by a Graphical user Interface which enhances the accessibility of the toolkit. The GUI framework is built using PyQt5, python binding of the cross-platform GUI toolkit Qt. The GUI framework allows users with minimal programming expertise to perform analysis.

  

#### 3.4. Terminology

  

#### 3.5. Launching the wrapper

iCOMIC can be launched using a simple command in the terminal.

```

$ ‘python mainwin_v32.py’

```

#### 3.6. Input file format

iCOMIC accepts input information in two different modes. The user can either feed the path to a folder containing raw fastq files or provide a table consolidating particulars of raw data.

  

>If you are uploading a folder of fastq files, all the files should be named in the specified format:

(sample_name)_(condition(tumor/normal))_Rep(replicate_number)_R(1 / 2).fastq

Example: hcc1395_normal_Rep1_R1.fastq

  

>If you choose upload from table mode, the sample information should be given in a tab delimited file with a header row.

The Column names should be:

`Sample` : The Sample name

`Unit` : The number of replicates

`Condition` : Nature of the sample, Normal or Tumor

`fq1` : The path of Read 1

`fq2` : Path of Read 2, if you are working with single-end reads only, the 'fq2' column can be left blank.

  

#### 3.7. Quick Guide

iCOMIC toolkit enables the ready analysis of RNA-Seq and Whole Genome/Exome Sequencing data. iCOMIC has an inbuilt library of tools with predefined valid combinations.

The user will have the freedom to choose any possible combination of tools. Figure 1 and 2 depicts the basic steps and outputs involved in DNA-Seq and RNA-Seq pipelines respectively.

![Figure 1: DNA Seq pipeline](https://github.com/anjanaanilkumar1289/iCOMIC_doc/blob/master/docs/dnaseq.png)


Figure 1: DNA Seq pipeline. This workflow indicates the analysis steps and major output files in DNA Seq pipeline.

![enter image description here](https://github.com/anjanaanilkumar1289/iCOMIC_doc/blob/master/docs/rnaseq.png)

Figure 2: The analysis steps and major output files in RNA Seq pipeline.

Here is a typical set of actions to run iCOMIC pipelines:
1. Select a pipeline.

2. Choose the mode of input

3. Input the required data fields.

4. Proceed to the next tab if you want to skip Quality Check.

5. Or click on the `Quality Control Results` button to view a consolidated MultiQC report of Quality statistics.

6. Check `yes` if you want to do trimming and also mention the additional parameters as per requirement.

7. Tool for Quality Control: FastQC

8. Tool for trimming the reads: Cutadapt

9. Choose the tools of interest from `Tool selection` tab and set the parameters as required

10. For the choice of aligner, the corresponding genome index file needs to be uploaded if available, or the user can generate the index file using the `Generate Index` button.

11. Click `Run` on the next tab to run the analysis.

12. If a warning button pops-up near the `Unlock` button, click on it to unlock the working directory.

13. Once the analysis is completed, `Results` tab will be opened.

14. DNA Seq results include a MultiQC report comprising the statistics of the entire analysis. A file consisting of the variants called and the corresponding annotated variant file.

15. Results for RNA Seq analysis include multiQC analysis statistics, R plots such as MA plot, Heatmap, PCA plot and box plot and list of differentially expressed genes.

#### 3.8. Output information

All outputs are stored in separate folders inside the main folder iCOMIC, for each pipeline along with log information.

- DNA-Seq
DNA Seq analysis generates five output folders as follows.
	- MultiQC
	Contains subfolders MultiQC, FastQC and Cutadapt. MultiQC contains consolidated `html` reports on the overall run statistics and a separate `html` file on merged FastQC reports of all the input samples. The folder FastQC contains quality reports of individual samples. It may also enclose FastQC report of trimmed reads if the user opts for trimming the input reads. The folder Cutadapt contains trimmed `fastq` files.
	- Aligner
Contain `bam` outputs generated by the aligner.
	- Variant Caller
This folder includes `vcf` files of identified variants.
	- Annotator
Contains annotated `vcf` files.
	- Index
This is an optional folder which contains the index files if the user chooses to generate index corresponding to the choice of aligner.
- RNA-Seq
RNA-Seq analysis generates five output folders inside the main folder iCOMIC, as follows.
	- MultiQC
Contains subfolders MultiQC, FastQC and Cutadapt. MultiQC contains consolidated `html` reports on the overall run statistics and a separate `html` file on merged FastQC reports of all the input samples. The folder FastQC contains quality reports of individual samples. It may also enclose FastQC report of trimmed reads if the user opts for trimming the input reads. The folder Cutadapt contains trimmed `fastq` files.
	- Aligner
Contain `bam` outputs generated by the aligner.
	- Expression Modeller
Contains count matrix representing the reads mapped to individual genes.
	- Differential expression
Contain a text file with a consolidated list of differentially expressed genes.
	- Index
This is an optional folder which contains the index files if the user chooses to generate index corresponding to the choice of aligner.
#### 3.9. FAQs

- What to do if run fails

- How to run iCOMIC?

- Dependences

- Installation issues
