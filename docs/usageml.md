## 13. Using iCOMIC for cTaG and NBDriver tools

#### 5.1. About cTaG

`cTaG` (classify TSG and OG) is a tool used to identify tumour suppressor genes (TSGs) and oncogenes (OGs) using somatic mutation data. The cTaG model returns the list of all genes labelled as TSG or OG or unlabelled along with predictions made by each model and whether the gene is among the top predictions.

##### 5.1.1. Input Requirements

A `maf` file is required to run the cTaG tool. A maf file can either be generated from the DNA-Seq output vcf file in the results tab or browsed locally. Added to that, you can mention the `parameters` required to run the cTag in the parameters option provided.

##### 5.1.2. Initialization of the analysis

You can click on the `run` button to initialize the analysis, once the necessary files have been uploaded.

##### 5.1.3. Results 

Once the analysis is completed, you can click on the `Results` button to view the results


#### 5.2. About NBDriver

`NBDriver` (NEIGHBORHOOD Driver) is a tool used to differentiate between driver and passenger mutations using features derived from the neighborhood sequences of somatic mutations. NBDriver returns a list of all mutations labelled as Driver or Passenger.

##### 5.2.1. Input Requirements

A `vcf` file is required to run NBDriver. The vcf file can either be browsed from the DNA-Seq output directory or locally.

##### 5.2.2. Initialization of the analysis

You can click on the `run` button to initialize the analysis, once the necessary files have been uploaded. 

##### 5.2.3. Results 

Once the analysis is completed, you can click on the `Results` button to view the results

