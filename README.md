# SUMMIT

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7034435.svg)](https://doi.org/10.5281/zenodo.7034435)

**Summary-level Unified Method for Modeling Integrated Transcriptome (SUMMIT)**, is a novel framework designed to improve the expression prediction model accuracy and the power of sequential TWAS by leveraging a very large expression quantitative trait loci (eQTL) summary-level dataset. One benefit of SUMMIT is that it can deal with genes with moderate-low expression heritability, which have been largely ignored by conventional TWAS.

### SUMMIT's workflow.

<p align="center">
<img src="illustration.bmp" alt="workflow" width="700"/>
</p>

Please cite the following manuscript for using SUMMIT and gene expression models built by SUMMIT:

> Zhang, Z., Bae, Y. E., Bradley, J. R., Wu, L., & Wu, C. (2021). SUMMIT: An integrative approach for better transcriptomic data imputation improves causal gene identification. medRxiv. Under Review.


### Replication:
For a complete replication of the results in our manuscript, please see our tutorials at osf.io (https://doi.org/10.17605/OSF.IO/BS3QU).

## Outline

1. [Training imputation model](#TRAIN)
2. [Association test](#TWAS)

## <a name="TRAIN"></a>Training models

```mainbody_cpp_pos_precise.R/mainbody_cpp_rsid_precise.R``` are SUMMIT's main functions to train imputation models. It has only 3 input arguments ```--name_batch``` ,```--method```, and ```wd```. ```--name_batch``` is the desired name for this batch of output and ```--method``` is designated penalized regression method that you can choose from (LASSO, ElNet, SCAD, MCP, MNet). ```--wd``` is the path of working directory.

### Data preparation

The Mainbody_cpp_xxx_precise functions require a very specific set of processed data to work properly. Unfortunately, we cannot share all the data with you due to privacy and confidentiality concern. Here, we provide a list of the data that we used. If you need further assistance on how the data were organized/processed, please reach out to us.

| Datasets referenced | Usage | How to obtain |
| ----- | ----- | ---- |
| gencode.v26.hg19.genes.rds | A look-up list for translating Ensembl IDs into gene names | downloaded from https://www.gencodegenes.org/human/release_26.html |
| Whole Blood_QCed_rpkm.rds | GTEx-7's subject's gene expression levels | downloaded from https://gtexportal.org and processed |
| response.8.RData | GTEx-8's subject's gene expression levels | downloaded from https://gtexportal.org and processed |
| subset-hapmap3 | eQTLGen summary statistics subsetted by HapMap3 | downloaded from https://www.eqtlgen.org/cis-eqtls.html with standard QC, then split into smaller files by gene name |
| seq.ref | Genotype matrix of reference panel (1000 Genomes Project) | downloaded from https://www.internationalgenome.org/data |
| seq.8 | Genotype matrix of GTEx-8 subjects | downloaded from https://gtexportal.org and process |
| chrX.OMNI.interpolated_genetic_map | Genetic distance of reference panel (1000 Genomes Project) | downloaded from https://github.com/joepickrell/1000-genomes-genetic-maps |

### Data alignment

SUMMIT offers two approaches to align the reference panel with the eQTL summary-level data. If you wish to align your data by rsID, use mainbody_cpp_rsid_precise.R. If you wish to align your data by chromosome+position, use mainbody_cpp_pos_precise.R.

### Example run

After we prepared the data, we can train imputation models via the following command.

```
Rscript mainbody_cpp_rsid.R \
--name_batch SCAD_1e-6_rsid \
--method SCAD \
--wd MY-WORK-DIR \
```

### Built-in parallel computing

Both ```mainbody_cpp_rsid_precise.R``` and ```mainbody_cpp_pos_precise.R``` contain a snippet that guarantees mutual exclusion for every subjob. Simply put, you can run ```mainbody_cpp_rsid_precise.R``` and ```mainbody_cpp_pos_precise.R``` using many parallel instances as you want and it will figure out if there is unfinished job on its own.

## <a name="TWAS"></a>Association test

### Gene expression prediction models built by SUMMIT

We have uploaded our pre-calculated expression imputation models (Tissue: whole blood) to osf.io (https://doi.org/10.17605/OSF.IO/7MXSA).

The osf.io repository contains two zip files. ``SUMMIT-weight-pos.zip`` contains models that use chromosome plus position (CHR + POS) to match the SNPs in our models to the SNPs in GWAS summary data; ``SUMMIT-weight-rsid.zip`` uses rsID to match.

### Pre-process your summary statistics using APSS

APSS is an interactive R function that can easily process GWAS summary statistics and shape GWAS summary statistics into any desired format.

APSS has 3 principal input arguments.

```directory.working``` is the working directory.

```filename``` is the name of the summary statistics file to be processed.

```BIG``` is the number of GBs and default is 2. If ```BIG``` is set as 2, then for any summary statistics file bigger than 2GB, APSS will do an exploratory read first. By doing so, APSS could significantly shorten the runtime and handle summary statistics files bigger than 10GB. 

### The must-have columns

You can use any summary statistics files with reasonable quality control just as long they contain specific columns for SUMMIT to work with.

The must-have columns are ```A1, A2, Z, CHR```.

If you are using the rsID-aligned models, you must also include ```SNP``` column.

If you are using the position-aligned models, you must also include ```POS``` column.

### The flags

```step2_pos.R``` is SUMMIT's function for association test using position-aligned models.

```step2_rsid.R``` is SUMMIT's function for association test using rsID-aligned models.

For both functions, the input arguments are:

```--models``` is the path of the pre-calculated models (Please make sure that the folder contains only model files).

```--path.ref``` is the path of the reference panel used plus the prefixes of reference panel files.

```--path.ss``` is the path of the reference panel (with "/" at the end).

```--trait``` is the path plus name of the input summary statistics files.

```--path.out``` is the path of the output folder (with "/" at the end).

```--parallel``` is the number of parallel instances

### Parallelization

Unlike other TWAS methods (e.g., PrediXcan), SUMMIT can be a bit more time-consuming. Because SUMMIT does not use a designated correlation matrix (LD matrix), the SUMMIT pipeline would spend more time matching models, summary statistics, and reference panel. In addition, each SUMMIT's pre-calculated model file contains models from 5 different types of penalized regression, hence more computation time is needed.

However, with proper parallelization, a complete association study usually takes less than 10 minutes. For example, if you decided to split the association tests into 20 smaller subjobs, an index (from 1 to 20) would need to be explicitly passed on to R from the global environment and the ```--parallel``` flag should be set to 20. Depending on your computing environment, you
may need to manually modify line 2 in ```step2_pos.R/step2_rsid.R```. 

### Example run

```
Rscript step2_pos.R \
--models my-model-folder-pos \
--path.ref my-reference-panel/1000G.EUR.ALLSNP.QC.CHR \
--path.ss my-ss-folder/ \
--trait my-trait-1 \
--path.out my-output-folder \
--parallel 100
```

```
Rscript step2_rsid.R \
--models my-model-folder-rsid \
--path.ref my-reference-panel/1000G.EUR.ALLSNP.QC.CHR \
--path.ss my-ss-folder/ \
--trait my-trait-2 \
--path.out my-output-folder \
--parallel 50
```

### Output format

| Column number   | Column name | Description |
| ----- | ----- | ---- |
| 1 | gene_symbol | Gene name  |
| 2 | gene_id | Ensembl ID  |
| 3 | chromosome | Chromosome |
| 4 | model_best | Best model |
| 5 | r2_test | Best model's $R^2$ on testing data  |
| 6-10 | p_ElNet | p-value of TWAS (Method is after the underscore) |
| 11-15 | z_ElNet | Z-score of TWAS (Method is after the underscore) |
| 16 | p_ACAT | The combined ACAT p-value |
| 17 | gene_pos | Gene position  |
| 18 | runtime | Runtime |

### Disclaimer

The built gene expression prediction models and software are provided “as is”, without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. in no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the models or the use or other dealings in the models.

### License

Maintainer: [Zichen Zhang] (zz17@fsu.edu)

[MIT](http://opensource.org/licenses/MIT)

Copyright (c) 2013-present, Zichen Zhang (zz17@fsu.edu), Chong Wu (cwu18@mdanderson.org)
