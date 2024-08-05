# RNAseqDEAnalysis
## Rscript for differential expression analysis of RNAseq data using DEseq2 after getting a merged gene count csv file

### If you do not have the merged_gene_counts.csv file, go to my preprocessing github to create it: https://github.com/tellmathieu/RNAseqPreprocess

Adapted from Patricia Baldrich's R script for differential expression analysis

## 1. Clone this repository to your working space

```
git clone https://github.com/tellmathieu/RNAseqDEAnalysis.git
```

## 2. Set up an environment to run this pipeline
- a) install mamba if you don't already have it - you'll have to go through the installation process or do a command like `module load mamba`
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
- b) Initialize mamba and set up an environment - this can take 5-10 minutes - I use specific versions of these programs to ensure everything works together. 
```
source ~/.bashrc #or you can log out and then log back in
mamba create -c conda-forge -c bioconda -y -n RNAseqDE r-upsetr=1.4.0 bioconductor-deseq2=1.42.0 bioconductor-apeglm=1.24.0 r-dplyr=1.1.4 r-ggplot2=3.5.1 r-ggrepel=0.9.5 r-rcolorbrewer=1.1_3 r-rlang=1.1.4 r-purrr=1.0.1 r-gplots=3.1.3.1 r-devtools=2.4.5
mamba activate RNAseqDE
```

## 3. Create metadata file (csv)

This is your description of all the different RNAseq samples you have, and the axes you want to compare when you run the differential expression analysis.

The two required columns are `sample` and `DE`. You can add more, but this program doesn't read them. The `sample` column must be the first column.

The `sample` column must match match the headers in your `merged_gene_counts.csv` file. The `DE` column should have your treatment names.
```
sample,DE
sample1,treatment1
sample2,treatment1
sample3,treatment1
sample4,treatment2
sample5,treatment2
sample6,treatment2
sample7,control
sample8,control
sample9,control
```

Here is an example of what my file looked like:
```
sample,treatment,replicate,DE
Col0_AWF_NT_R4,Col0_AWF_NT,R4,Col0_AWF_NT
Col0_AWF_NT_R5,Col0_AWF_NT,R5,Col0_AWF_NT
Col0_AWF_NT_R6,Col0_AWF_NT,R6,Col0_AWF_NT
Col0_CL_R1,Col0_CL,R1,Col0_CL
Col0_CL_R2,Col0_CL,R2,Col0_CL
Col0_CL_R3,Col0_CL,R3,Col0_CL
Col0_CL_ribo_R1,Col0_CL_ribo,R1,Col0_CL_ribo
Col0_CL_ribo_R2,Col0_CL_ribo,R2,Col0_CL_ribo
Col0_CL_ribo_R3,Col0_CL_ribo,R3,Col0_CL_ribo
```

## 4. Create treatment vs comparisons (treatmentComparisons) file (csv)

The `num` column is a sequentially numbered column.

The `treatment` column are the treatments that you treated your samples with. These must be exactly the same as the `DE` column in step 3.

The `control` column are what you consider controls in your experiment. These must be exactly the same as the `DE` column in step 3.

```
num,treatment,control
1,treatment1,control
2,treatment2,control
```

Here is an example of my file that I used:
```
num,treatment,control
1,Col0_AWF_NT,Col0_CL
2,Col0_LSW,Col0_CL
3,Col0_LSW,Col0_AWF_NT
4,Col0_CL_ribo,Col0_CL
```

## 5. Check that you have your input files and input variables

At this point, you should have 3 `.csv` files:
- `merged_gene_counts.csv`
- `sampleMetadata.csv`
- `treatmentComparisons.csv`

Additionally, you need 3 more variables to run this Rscript:
- alpha: The alpha value for significance in the statistical analysis. Commonly, "0.05" is a good starting point, and you can adjust if you get too few or too many differentially expressed genes.
- resultsPath: The filepath to your results folder. This folder must exist before you run the Rscript. An example would be "home/DEseq/results."
- ctrl: The name of the control treatment (must exactly match a treatment from the `DE` column in step 3). This is the overall control for the experiment, which is good to have when running DEseq2


## 6. Run Rscript

Fill in the items with `{}` brackets.
```
Rscript deseq.R \
	{merged_gene_counts.csv} \
	{sampleMetaData.csv} \
	{alpha} \
	{resultsPath} \
	{treatmentComparisons.csv} \
	{ctrl}
```

Here is an example of an actual command:
```
Rscript deseq.R \
	example/merged_gene_count.csv \
	example/sampleMetaData.csv \
	0.05 \
	home/tmathieu/RNAseqDEAnalysis/results \
	example/treatmentComparisons.csv \
	Col0_CL
```




