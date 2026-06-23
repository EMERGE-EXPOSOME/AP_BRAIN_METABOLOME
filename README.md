# AP_BRAIN_METABOLOME

Code respiratory for the "Air pollution brain metabolomics study" led by Dr. Mengyuan Ren.

## Brief description about this repo
This repository contains code files, and source data. This study, led by Dr. Ren et al, aims to explore air pollution associated brain metabolic perturbtaions in the context of Alzheimer's disease. This study investigates the associations between long-term air pollution exposure and brain metabolome, using metabolome-wide association study (MWAS), pathway enrichment analysis to uncover potential biological mechanisms underlying these associations, and metabolite annotation to confirm metabolites of Level-1 evidence. Using mixture analysis, this study uncovers joint effects of air pollution on brain metabolites. 

## 1. Equipments that used for running code files
### 1.1 Hardware
Below is the hardware that used for the manuscript. We recommend a hardware with at least 4GB RAM and 4 cores.<br/>

Dell Latitude 7350, Intel(R) Core(TM) Ultra 7 165U 1.70GHz <br/>
RAM: 32GB <br/>
System type: Windows 11 (64-bit) operating system <br/>

### 1.2 Software (versions)
- R (version 4.4.1, open-source, GNU General Public License) <br/><br/>

*The code files were also independently tested on another device with parameters below:* <br/>

#### Hardware:
Dell laptop, Intel(R) Core(TM) Ultra 5 135U (2.10 GHz) <br/>
RAM: 16.0 GB <br/>
System type: Windows 11 Enterprise, 64-bit operating system, x64-based processor <br/>

#### Software:
- R version 4.5.0 (2025-04-11 ucrt) <br/><br/>


For R package versions, see below: 
| Package name | Package version |
| ------------- | ------------- |
| apLCMS | 6.6.9
| ComplexUpset | 1.3.3
| corrplot | 0.95
| data.table | 1.16.4
| ddpcr | 1.15.2 
| dplyr | 1.1.4
| ggbreak | 0.1.4
| ggplot2 | 3.5.1
| ggpubr | 0.6.0
| ggridges | 0.5.6
| ggtext | 0.1.2
| glue | 1.8.0
| kableExtra | 1.4.0
| lubridate | 1.9.4
| magrittr | 2.0.3
| metapone | 1.12.0
| patchwork | 1.3.0
| purrr | 1.0.4
| qgcomp | 2.16.4
| R6 | 2.6.1
| RColorBrewer | 1.1-3
| readxl | 1.4.3
| rlist | 0.4.6.2
| sf | 1.0-19
| statmod | 1.5.0
| stringr | 1.5.1
| tableone | 0.13.2
| tictoc | 1.2.1
| tidyverse | 2.0.0
| UpSetR | 1.4.0
| xMSanalyzer | 2.0.6.1

## 2. Instructions for R package installation
### 2.1 Install R package pacman and BiocManager, both of which can be used to install multiple packages in one shot
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("pacman", quietly = TRUE))
    install.packages("pacman")
```


### 2.2 Install packages indicated in code files, taking "tidyverse" and "ggplot2" as example: 
```r
pacman::p_load("tidyverse", "ggplot2")
```


In some cases, installing package may get failed for several reasons. Here we listed instructions for installing specific packages, including xMSanalyzer, and apLCMS, both of which were used for metabolite annotation in this study. 


Before installing these two R packages, you may need to install ROCS package using .tar.gz file ([see here](https://github.com/EMERGE-EXPOSOME/OMIT_PTB_AIM2/blob/main/ROCS_1.3.tar.gz)), as well as other dependencies, including but are not limited to: xcms, ProGenerics, sva, mixOmics, WGCNA, GO.db


```r
BiocManager::install(c("xcms", "ProGenerics", "sva", "mixOmics", "WGCNA", "GO.db"))
```


When you finish above, run the remote::install code below finishing installation
```r
remotes::install_github("kuppal2/xMSanalyzer")
remotes::install_github("tianwei-yu/apLCMS")
```


It usually takes **2-3 hours** installing all R packages as well as their dependencies. 


## 3. Run code on demo dataset

xxxxxxxxxxxxx


## 4. How to run the software on demo data

There are several integrated development environment (IDE) widely used to run R code, including Rstudio, vscode, etc. To run R code on Rstudio, please [see here](https://docs.posit.co/ide/user/ide/guide/code/execution.html) for details. To run R code on vscode, please [see here](https://code.visualstudio.com/docs/languages/r) for instructions. We **recommend using Rstudio** owing to its user-friendly interface. 


## Contact

If you have any questions, please reach out to Drs. Mengyuan Ren (mren30@emory.edu), Anke Huels, and Donghai Liang. 









