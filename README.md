# Reliability analysis of a ship fender using a vine-copula

## Create R Environment

Instructions for running analysis with local Anaconda installation. The following steps should be run in a terminal and will install a conda environment for R and run the analysis.


```
conda create --name r-vine-fender
conda activate r-vine-vine-fender
conda install -c conda-forge r-base=4.1.3
conda install --file requirements.txt
Rscript packs.r
Rscript run_reliability.r
```
