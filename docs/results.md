## 8. Viewing and analyzing results

#### 8.1. MultiQC reports
The tool multiQC compiles the analysis statistics for different tools and provides a consolidated report. The tool is used to visualize the analysis results at multiple stages in iCOMIC. The Quality Control part in iCOMIC analyses the quality of all input reads using FastQC. MultiQC compiles the FastQC reports for each sample and provides a consolidated comprehensive report. MultiQC is used again to summarise the entire analysis statistics. In the case of Whole genome sequencing, the MultiQC report includes a compiled FastQC report, alignment statistics and statistics of the variants identified. On the other hand, the MultiQC report in the RNA-Seq analysis part includes results from tools such as FastQC, Cutadapt, and STAR.

#### 8.2. Plots generated in RNA-Seq
Differential Expression tools generate R plots such as MA plot, Heatmap, PCA plot and box plot displays the predicted differentially expressed genes. MA plot helps to find log2 fold changes,Heatmap helps in exploring the count matrix, PCA Plot visualizes the overall effect of experimental covariates and batch effects and box plots used to find count outliers.
#### 8.3. List of differentially expressed genes in RNA-Seq
#### 8.4. Variants Called in DNA-Seq
iCOMIC displays the variants identified in `vcf` format. In the results tab, the user can click on the `click` button and a pop-up with the `vcf` file will be displayed. 
#### 8.5. Annotated variants in DNA-Seq
Here the `vcf` file of annotated variants are displayed.
