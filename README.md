R Environment

for vine copula example

[This](https://andrewpwheeler.com/2022/04/08/managing-r-environments-using-conda/) is a useful post.

Find newest available r package on conda, like [`r-base`](https://anaconda.org/conda-forge/r-base). There are several options...but be sure to check the version for your OS! Seems like windows lags behind the others, especially Linux. I chose 4.1.3 on Feb 14, 2023.
```
conda create --name r-vine-fender
conda activate r-vine-vine-fender
conda install -c conda-forge r-base=4.1.3
conda install --file requirements.txt
Rscript packs.r
Rscript run_reliability.r
```


Had to extract the packages that are in the base install, then installed from CRAN using `packs.R`. Everything seemed to be ok. Loading everything with `libraries.r` ran without error, though.

