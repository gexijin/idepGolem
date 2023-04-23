# iDEP (integrated Differential Expression & Pathway analysis)

## Description

iDEP is an bioinformatics tool for analyzing gene expression data graphically, interactively and reproducibly. Hosted at [http://bioinformatics.sdstate.edu](<http:://bioinformatics.sdstate.edu/idep11/>), iDEP is developed as an R package based on the [Golem framework](https://thinkr-open.github.io/golem/).

## Run iDEP locally as an R package or via a Docker container
After redesigning our database, the new iDEP can be easily run on your laptop or a local server. Please note that iDEP is being updated frequently. If you run locally, please upgrade to our most recent version on a monthly basis. 
After being installed as an R package, iDEP can be started from R using the ```iDEP::run_app()``` command. This method is time-consuming as it requires the installation of all of the 355 dependent R packages, which takes about an hour on Windows 10. It can take much more time for Mac or Linux users, as many of these R packages require trouble-shooting. 
A faster alterative is to use the [iDEP Docker image](https://hub.docker.com/repository/docker/gexijin/idep/general) on DockerHub. The only requirement is install Docker software on Linux, and Docker Desktop on Windows or MacOS.

### System requirements
Most of modern laptop can run iDEP locally. Minimum storage 10GB. Minimum Memory 4GB. 

### Windows: Docker Desktop (~20 minutes, recommended)
1. Install [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)
2. Start Windows PowerShell as an **Administrator**. From Windows search bar, type **PowerShell** to find the Windows PowerShell app. And then select **Run as an Administrator**. See [here](https://www.howtogeek.com/742916/how-to-open-windows-powershell-as-an-admin-in-windows-10/) for detailed steps. 
3. [Enable WSL2](https://learn.microsoft.com/en-us/windows/wsl/install-manual)
4. Start the Docker app. From Windows search bar, type **Docker**, and then select Run. Accept the terms when asked. The Docker engine is now running in the background.
5. Pull the iDEP Docker image from PowerShell command line.
```console
docker pull gexijin/idep:latest
```
5. Run docker container from PowerShell.
```console
docker run --pull -d --name idep -p 3838:3838 gexijin/idep:latest 
```
6. You can now use iDEP locally by starting your web browser and enter localhost:3838 in the web address.

Note that the Docker container needs to be kept running as a webserver. If you restart your computer or accidentally closed the Command Prompt window, you need to re-run Step 5. 

### Windows: installation as an R package (~1 hour)
1. Install a [recent version of R](https://cloud.r-project.org/) such as R 4.30. 
2. Optional: Install IDE such as [RStudio](https://posit.co/download/rstudio-desktop/) or [VS Code](https://code.visualstudio.com/).
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

### Windows: Developer mode (~1 hour)
This method enables users to customize the iDEP by changing the source code. We also welcome contributions by submitting Pull Requests.
1. Install a [recent version of R](https://cloud.r-project.org/) such as R 4.30. 
2. Optional: Install IDE such as [RStudio](https://posit.co/download/rstudio-desktop/) or [VS Code](https://code.visualstudio.com/).
3. Obtain a copy of the source code. This can be done manually by downloading an zip file by clicking on the green Code button above. Alternatively, you can fork this repository, install GitHub Desktop, and clone this repository locally using an URL(https://github.com/gexijin/idepGolem.git). 
4. Start R and install the golem package.
```{R}
install.packages("golem")
```
5. Open the downloaded reporsitory from R. The Shiny app can be started by running the ```run_dev.R``` script. If R packages are missing, install them and try again. The app is devided into 11 Shiny Modules.


### MacOS: Docker container (~10 minutes on MacBook Air)
1. Download Docker Desktop follow instructions [here](https://www.docker.com/products/docker-desktop/). Make sure you choose the correct version based on your CPU type. For MacBook Air, I chose ```Apple Chip```. If you use an Mac computer that uses Intel CPU, choose the ```Intel Chip```.
2. Install Docker Desktop. First double-click the downloaded Docker.dmg file in the Downloads folder. Drag the Docker icon into the Application folder from the pop-up window. From Lunch pad, or the Application folder, click on the Docker icon. Click ```Open``` when asked. Accept the Terms and the Docker engine is running.
3. Start a Terminal window by clicking the Launchpad, and type ```terminal``` in the search field. Then click the Terminal app. 
4. Pull the iDEP Docker image using this command in the Terminal window.
```console
docker pull gexijin/idep:latest
```
5. Run docker container from Linux command line.
```console
docker run --pull -d --name idep -p 3838:3838 gexijin/idep:latest 
```
6. You can now use iDEP locally by starting your web browser and enter localhost:3838 in the web address.

Note that the Docker container needs to be kept running as a webserver. If you restart your computer or accidentally closed the Command Prompt window, you need to re-run Step 5. 


### Linux: Docker container (~10 minutes)
1. Install Docker Engine follow instructions [here](https://docs.docker.com/engine/install/). Many Linux systems have Docker installed by default. On Ubuntu, I used the following scripts.
```console
cd
curl -fsSL https://get.docker.com -o install_docker.sh
sudo sh install_docker.sh
```
2. Pull the iDEP Docker image from the Linux command line.
```console
docker pull gexijin/idep:latest
```
3. Run docker container from Linux command line.
```console
docker run --pull -d --name idep -p 3838:3838 gexijin/idep:latest 
```
4. You can now use iDEP locally by starting your web browser and enter localhost:3838 in the web address. If the is an internal server, replace the localhost with the IP address and make sure that the port 3838 is open.


