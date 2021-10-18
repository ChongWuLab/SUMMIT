# SUMMIT

SUMMIT (Summary-level Unified Method for Modeling Integrated Transcriptome).

More details of SUMMIT is described in,

> [LINK-1]

## Outline

1. [Training imputation model](#TRAIN)
2. [Association test](#TWAS)

## <a name="TRAIN"></a>Training models

mainbody_cpp_final.R is SUMMIT's main function to train imputation models. It has only 2 input arguments ```--name_batch``` and ```--method```. ```--name_batch``` is the desired name for this batch of output and ```--method``` is designated penalized regression method that you can choose from (LASSO, ElNet, SCAD, MCP, MNet).

### Data preparation

#### Reference panel

We used data from 1000 Genome European ancestry as our reference panel. With standard quality control procedure, we splitted the processed data by chromosome.

#### eQTLgen summary statistics

As for the summary-level eQTLgen data, you can download from https://www.eqtlgen.org/cis-eqtls.html. Notice that the eQTLgen dataset is over 3.5Gb, we have splitted them into small ```.RData``` files (You can download from LINK-2) by gene_id after standard quality control. If you do wish to process the summary-level eQTLgen data on your own, please do not change the column names or column order from eQTL-gen's original data.

#### GTEx-7 and GTEx-8 data

SUMMIT used 


### Alignment

SUMMIT offers two approaches to align the reference panel with the summary-level data.

### Example run

After we prepared the data, we can run SUMMIT via the following single line.

```
Rscript mainbody_cpp_rsid_precise.R \
--name_batch SCAD_1e-6_rsid \
--method SCAD \
```

### Trained imputaion models

All the ready-to-use SUMMIT Whole_Blood imputation models can be found here: LINK-5.

## <a name="TWAS"></a>TWAS

### Pre-process your summary statistics

APSS is a interactive R function to easily process GWAS summary statistics and shape GWAS summary statistics into any desired format. 

### The must-have columns

The must-have columns for SUMMIT are A1, A2, Z, CHR. If you are using the rsid-aligned models, you must include SNP column. If you are using the position-aligned models, you must include POS column.

### Example run

```
Rscript step2_pos.R \
--path.ref my-reference-panel/1000G.EUR.ALLSNP.QC.CHR \
--trait my-trait-1 \
--path.out my-output-folder \
```

```
Rscript step2_rsid.R \
--path.ref my-reference-panel/1000G.EUR.ALLSNP.QC.CHR \
--trait my-trait-1 \
--path.out my-output-folder \
```

## License

Maintainer: [Zichen Zhang] (zz17@my.fsu.edu)

[MIT](http://opensource.org/licenses/MIT)

Copyright (c) 2013-present, Zichen Zhang (zz17@my.fsu.edu)
