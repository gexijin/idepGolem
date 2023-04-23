# iDEP (integrated Differential Expression & Pathway analysis)

## Description

This is the Shiny app of [iDEP](<http:://bioinformatics.sdstate.edu/idep11/>), developed as an R package based on the Golem framework.
iDEP is a online analytics tool that enables uers to analyze gene expression data graphically, interactively and reproducibly. 


## Installation

### Windows. Method 1. installation of R package
This takes about one hour. See below for running on Windows using Docker Desktop.
1. Install a [recent version of R](https://cloud.r-project.org/) such as R 4.30. 
2. Optional: Install [RStudio](https://posit.co/download/rstudio-desktop/) 
3. Start R and run these commands. It takes about an hour to install the 355 dependencies! If there are issues with the installation of some of the individual packages, you need to resolve the issues and try to install them.

```{R}
install.packages("devtools")
devtools::install_github("https://github.com/gexijin/idepGolem")
```

You might get the following warnings, which could be ignored.

```{R}
WARNING: Rtools is required to build R packages, but is not currently installed.
package ‘KEGG.db’ is not available for this version of R
A version of this package for your version of R might be available elsewhere,
2: In i.p(...) : installation of package ‘GO.db’ had non-zero exit status
3: In i.p(...) :
  installation of package ‘/tmp/RtmpCpiUsZ/file2326c2a50656f/PGSEA_1.60.0.tar.gz’ had non-zero exit status
```
4. Start the iDEP from R
```{R}
idepGolem::run_app()
```

### Windows, using Docker Desktop
1. Install [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)
2. [Enable WSL2](https://learn.microsoft.com/en-us/windows/wsl/install-manual)
3. Start a Command Prompt. Press the Win + R keys on your keyboard, then type cmd, and then click OK.
4. Pull the iDEP Docker image from the Command Prompt
```
docker pull gexijin/idep:latest
```
5. Run docker container from Command Prompt
```
docker run --pull -d --name idep -p 3838:3838 gexijin/idep:latest 
```
6. You can now use iDEP locally by starting your web browser and enter localhost:3838 in the web address.

Note that step 5 needs to be repeated if you restart your computer or closed the Command Prompt window.

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
