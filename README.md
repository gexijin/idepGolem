# iDEP (integrated Differential Expression & Pathway analysis)

## Description

iDEP is an bioinformatics tool for analyzing gene expression data graphically, interactively and reproducibly. Hosted at [http://bioinformatics.sdstate.edu](<http:://bioinformatics.sdstate.edu/idep11/>), iDEP is developed as an R package based on the Golem framework.

## Run iDEP locally by installing the R package or using a Docker container
After redesigning our database, the new iDEP can be easily run on your laptop or a local server. Please note that iDEP is being updated frequently. If you run locally, please upgrade to our most recent version on a monthly basis. 
After being installed as an R package, iDEP can be started from R using the ```iDEP::run_app()``` command. This method is time-consuming as it requires the installation of all of the 355 dependent R packages, which takes about an hour on Windows 10. It can take much more time for Mac or Linux users, as many of these R packages require trouble-shooting. 
A faster alterative is to use the [iDEP Docker image](https://hub.docker.com/repository/docker/gexijin/idep/general) on DockerHub. The only requirement is install Docker software on Linux, and Docker Desktop on Windows or MacOS.

### Windows installation as an R package
1. Install a [recent version of R](https://cloud.r-project.org/) such as R 4.30. 
2. Optional: Install [RStudio](https://posit.co/download/rstudio-desktop/) 
3. Start R and run these commands. It takes about an hour to install the 355 dependencies! If there are issues with the installation of some of the individual packages, you need to resolve the issues and try to install them.

```{R}
install.packages("devtools")
devtools::install_github("https://github.com/gexijin/idepGolem")
```
You might get the following warnings, which can be ignored.

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

### Windows using Docker Desktop
1. Install [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)
2. [Enable WSL2](https://learn.microsoft.com/en-us/windows/wsl/install-manual)
3. Start a Command Prompt on Windows. Press the Windows key + R keys on your keyboard. Type cmd, and then click OK. For other methods see [here](https://www.howtogeek.com/235101/10-ways-to-open-the-command-prompt-in-windows-10/).
4. Pull the iDEP Docker image from the Command Prompt
```
docker pull gexijin/idep:latest
```
5. Run docker container from Command Prompt
```
docker run --pull -d --name idep -p 3838:3838 gexijin/idep:latest 
```
6. You can now use iDEP locally by starting your web browser and enter localhost:3838 in the web address.

Note that the Docker container needs to be kept running as a webserver. If you restart your computer or accidentally closed the Command Prompt window, you need to re-run Step 5. 

### Linux using Docker
1. Install Docker Desktop follow instructions [here](https://docs.docker.com/engine/install/). Many Linux systems have Docker installed by default. On Ubuntu, I used the following scripts.
```
cd
curl -fsSL https://get.docker.com -o install_docker.sh
sudo sh install_docker.sh
```
2. Pull the iDEP Docker image from the Command Prompt
```
docker pull gexijin/idep:latest
```
3. Run docker container from Command Prompt
```
docker run --pull -d --name idep -p 3838:3838 gexijin/idep:latest 
```
6. You can now use iDEP locally by starting your web browser and enter localhost:3838 in the web address. If the is an internal server, replace the localhost with the IP address and make sure that the port 3838 is open.


