# iDEP (integrated Differential Expression & Pathway analysis)

## Description

This is the Shiny app of [iDEP](<http:://bioinformatics.sdstate.edu/idep11/>), developed as an R package based on the Golem framework.
iDEP is a online analytics tool that enables uers to analyze gene expression data graphically, interactively and reproducibly. 


## Installation

### Windows (tested on Windows 10)
1. Install a [recent version of R](https://cloud.r-project.org/) such as R 4.30. 
2. Optional: Install [RStudio](https://posit.co/download/rstudio-desktop/) 
3. Start R and run these commands. It takes about an hour to install the 355 dependencies!

```{R}
install.packages("devtools")
devtools::install_github("https://github.com/gexijin/idepGolem")
```

You might get the following warnings, which could be ignored.

```{R}
WARNING: Rtools is required to build R packages, but is not currently installed.

* DONE (idepGolem)
Warning messages:
1: package ‘KEGG.db’ is not available for this version of R

A version of this package for your version of R might be available elsewhere,
see the ideas at
https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages 
2: In i.p(...) : installation of package ‘GO.db’ had non-zero exit status
3: In i.p(...) :
  installation of package ‘/tmp/RtmpCpiUsZ/file2326c2a50656f/PGSEA_1.60.0.tar.gz’ had non-zero exit status
```

### Linux (tested with Ubuntu v. 22)
1. Install R version 4.30  following [instructions](https://cloud.r-project.org/bin/linux/ubuntu/). 
2. Install [cMake](https://cmake.org/install/), required by the R package nloptr.
```
wget https://github.com/Kitware/CMake/releases/download/v3.26.3/cmake-3.26.3.tar.gz
tar -xzvf cmake-3.26.3.tar.gz
cd cmake-3.26.3
./bootstrap
make 
make install
```
3. Install [libproj-dev](), required by the ggalt package.
```
sudo apt update
sudo apt install libproj-dev
```
4. Install other Linux packages
```
sudo apt install -y \
libcurl4-openssl-dev \
libxml2-dev \
libxml2  \
libssl-dev \
libudunits2-dev \
libmariadbclient-dev \
libpng-dev
```
5. Start R and run these:
```{R}
install.packages("devtools")
devtools::install_github(
  "https://github.com/gexijin/idepGolem/tree/data107",
  upgrade =  "never"
)
```


## Run iDEP locally

From the R console:
```{R}
> idepGolem::run_app() 
Loading required package: shiny

Listening on http://127.0.0.1:7045
Browsing http://127.0.0.1:7045
```
