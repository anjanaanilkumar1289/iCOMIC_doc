## 10. Tutorials

#### 10.1. What is it?
`iCOMIC` (Integrating Context Of Mutation In Cancer) is an open-source, standalone tool for genomic data analysis characterized by a Python-based Graphical User Interface, automated bioinformatics pipelines for analyzing `Whole genome/exome`  and `transcriptomic` data along with Machine Learning tools, `cTaG` and `NBDriver` for cancer related data analysis. It serves as a point and click application facilitating genomic data analysis accessible to researchers with minimal programming expertise. iCOMIC takes in raw sequencing data in FASTQ format as an input, and outputs insightful statistics on the nature of the data. iCOMIC toolkit is embedded in `‘Snakemake’`, a workflow management system and is characterized by a user-friendly GUI built using `PyQt5` which improves its ease of access. The toolkit features many independent core workflows in both whole genomic and transcriptomic data analysis pipelines. nes.

#### 10.2. Prerequisites:
- Linux/Windows/Mac platform
- Python 3.6 and above
- Miniconda
- iCOMIC package downloaded from GitHub
- Memory requirement

#### 10.3. Installation

Installation is easy as we provide a requirements.txt file comprising all the software dependencies. Once you clone the iCOMIC github repository, you can install all the associated dependencies using the command below. Every additional software requirement will be managed by the conda environment.

`$ pip install requirements.txt`

#### 10.4. Testing

#### 10.5. Analysis quick guide:

##### 10.5.1 Launching the wrapper

iCOMIC can be launched using a simple command in the terminal.

`$ ‘python iCOMIC_v01.py’`

##### 10.5.2 Running iCOMIC: A quick walkthrough

Here is a typical set of actions to run iCOMIC pipelines:
- Select a pipeline
- Choose the mode of input
- Input the required data fields
- Proceed to the next tab if you want to skip Quality Check
- Or click on the Quality Control Results button to view a consolidated MultiQC report of Quality statistics
- Check yes if you want to do trimming and also mention the additional parameters as per requirement
- Tool for Quality Control: FastQC
- Tool for trimming the reads: Cutadapt
- Choose the tools of interest from Tool selection tab and set the parameters as required
- For the choice of aligner, the corresponding genome index file needs to be uploaded if available, or the user can generate the index file using the Generate Index button
- Click Run on the next tab to run the analysis
- Once the analysis is completed, the Results tab will be opened
- DNA-Seq results include a MultiQC report comprising the statistics of the entire analysis, a file consisting of the variants called and the corresponding annotated variant file
- Results for RNA-Seq analysis include multiQC analysis statistics, R plots such as MA plot, Heatmap, PCA plot and box plot and list of differentially expressed genes
- Proceed to cTaG/NBDriver tab for further analysis if needed

##### 10.5.3 Adding samples: step one

iCOMIC accepts input information in `two` different modes. In the `first` method, the user can feed the path to a folder containing raw fastq files. For the direct upload of a sample folder, the folder should contain only the samples and the sample file names should be in a specific format:

`{sample_name}_{condition}_Rep{replicate_number}_R{1 / 2}.fastq`

- `{sample_name}` should be replaced with the sample name
- `{condition}` should be replaced with the nature of the sample, normal or tumor. If you are using a germline variant calling pipeline, the condition should be normal for all the samples
- `{replicate_number}` should be replaced by the number of replicate
- If the sample is paired end, `{1 / 2}` should be replaced by 1 or 2 accordingly for forward and reverse sequences. If the sample is single end, {1 / 2} can be replaced by 1


##### 10.5.4 Adding samples: step two
##### 10.5.5 Adding samples: specifying DNA-seq workflow
##### 10.5.6 Adding samples: specifying RNA-seq workflow
##### 10.5.7 Review of Sample quality
##### 10.5.8 Specifying analysis settings DNA seq
##### 10.5.9 Setting up differential gene expression analysis
##### 10.5.10 Submitting the analysis
##### 10.5.11 Retrieving the data
##### 10.5.12 Analysis with BAM input
##### 10.5.13 Running cTaG
##### 10.5.14 Running NBDriver

#### 10.6 Retrieving logs:
