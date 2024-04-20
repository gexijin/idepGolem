## Integrated Differential Expression & Pathway analysis (iDEP)
[Original repository](https://github.com/iDEP-SDSU/idep)
## Description

[iDEP](http://bioinformatics.sdstate.edu/idep11/)  is a bioinformatics platform for analyzing gene expression data graphically, interactively, and reproducibly. The input file is a gene-level expression matrix derived from RNA-Seq, microarray, proteomics, or other methods. Hosted at [http://bioinformatics.sdstate.edu](<http:://bioinformatics.sdstate.edu/idep11/>), iDEP is developed as an R package based on the [Golem framework](https://thinkr-open.github.io/golem/), by a small team led by Dr. [Steven Ge](https://twitter.com/StevenXGe). See [documentation](https://idepsite.wordpress.com/) and [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6). 

## License
(CC BY-NC 3.0) Non-commercial use. For local installation at private institutions, please [contact us](https://www.sdstate.edu/directory/xijin-ge).

## Run iDEP locally on [Windows](https://github.com/gexijin/idepGolem#windows-docker-desktop-20-minutes-recommended), [MacOS](https://github.com/gexijin/idepGolem#macos-docker-container-10-minutes-on-macbook-air), or [Linux](https://github.com/gexijin/idepGolem#linux-docker-container-10-minutes)
Following the redesign of our database, the updated iDEP can now be effortlessly executed on a laptop or local server. Please be aware that iDEP undergoes frequent updates; if you are running the software locally, we recommend updating to the most recent version on a monthly basis.

Once installed as an R package, iDEP can be initiated from R using the ```iDEP::run_app()``` command. However, this approach can be time-consuming, as it necessitates the installation of all 355 dependent R packages. This process takes approximately an hour on Windows 10, and potentially even longer for Mac or Linux users due to potential troubleshooting requirements for many of the R packages.

A more efficient alternative is to utilize the [iDEP Docker image](https://hub.docker.com/repository/docker/gexijin/idep/general) available on DockerHub. The only prerequisite is to install the Docker software on Linux or Docker Desktop on Windows or MacOS.

Please note that local installation is free for non-profit organizations only.

### System requirements
Most of modern laptop can run iDEP locally, if it has more than 10GB storage and 4GB of memory. 

### Windows: Docker Desktop (~20 minutes, recommended)
Just follow this detailed [video](https://youtu.be/EJiNG9uUq5g). No prior experience with Docker is needed. 
1. Install [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/).
2. Start Windows PowerShell as an **Administrator**. From Windows search bar, type **PowerShell** to find the Windows PowerShell app. And then select **Run as an Administrator**. For details, see [here](https://www.howtogeek.com/742916/how-to-open-windows-powershell-as-an-admin-in-windows-10/). 
3. Enable [Windows Subsystem for Linux 2 (WSL2).](https://learn.microsoft.com/en-us/windows/wsl/install-manual)
4. Start the Docker app. From Windows search bar, type **Docker**, and then select **Run as an Administrator**. Accept the terms when asked. The Docker engine is now running in the background.
5. Pull the iDEP Docker image from DockerHub and start a container from PowerShell.

```console
docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest 
```
6. You can now use iDEP locally by starting your web browser and enter **localhost:3838** in the address bar.

Note that the Docker engine is now running in the backgroup, acting as a webserver. It works even if you restart your computer. To stop it, run these from Windows PowerShell: 
```console
docker stop idep 
docker rm idep
```
After stopping it, you can restart it by repeating Step 5, which also pulls the latest Docker image from DockerHub. Make sure you update your image at least on a monthly basis.

### Windows: iDEP as an R package (~1 hour)
1. Install a [recent version of R](https://cloud.r-project.org/). If your R is a few years old, uninstall it, and manually delete the folder that contains all existing R packages.
2. Optional: Install IDE such as [RStudio](https://posit.co/download/rstudio-desktop/) or [VS Code](https://code.visualstudio.com/).
3. Start R and run these commands. It takes about an hour to install the 355 dependencies! If there are issues with the installation of some of the individual packages, you need to resolve the issues and try to install them.

```{R}
install.packages("devtools")
devtools::install_github("https://github.com/gexijin/idepGolem", upgrade = "never")
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
4. Start iDEP from R command console.
```{R}
idepGolem::run_app()
```
The benefit of this approach is that you always have the most recent version from GitHub. The next time you install iDEP, it will take much less time as the dependencies have been installed. 

### Windows: Developer mode (~1 hour)
With this method, users can customize iDEP by changing the source code. We also welcome contributions through Pull Requests on GitHub.
1. Install a [recent version of R](https://cloud.r-project.org/), such as R 4.30. 
2. Optional: Install IDE such as [RStudio](https://posit.co/download/rstudio-desktop/) or [VS Code](https://code.visualstudio.com/).
3. Obtain a copy of the source code. This can be done manually by downloading an zip file by clicking on the green Code button above. Alternatively, you can fork this repository, install GitHub Desktop, and clone this repository locally using an URL(https://github.com/gexijin/idepGolem.git). 
4. Start R and install the golem package.
```{R}
install.packages("golem")
```
5. Optional: Install all dependencies. We pretend to install the idepGolem package. This take about an hour.
```{R}
install.packages("devtools")
devtools::install_github("https://github.com/gexijin/idepGolem", upgrade = "never")
remove.packages("idepGolem")
```
6. Open the downloaded reporsitory from R. The Shiny app can be started by running the ```run_dev.R``` script. If some R packages are missing, install them and try again. Note that the app is devided into 11 Shiny Modules.

### MacOS: Docker container (~10 minutes on MacBook Air)
See [video](https://youtu.be/u8Gdog4VAGc) for more details. No prior experience with Docker is needed.
1. Download Docker Desktop follow instructions [here](https://www.docker.com/products/docker-desktop/). Make sure you choose the correct version based on your CPU type. For MacBook Air, I chose **Apple Chip**. If you use an Mac computer with Intel CPU, choose the **Intel Chip** instead.
2. Install Docker Desktop. First double-click the downloaded Docker.dmg file in the Downloads folder. Drag the Docker icon into the Application folder from the pop-up window. From Lunch pad, or the Application folder, click on the Docker icon. Click **Open** when asked. Accept the Terms and the Docker engine is running.
3. Start a Terminal window by clicking the Launchpad, and type ```terminal``` in the search field. Then click the Terminal app. 
4. Pull the iDEP Docker image from DockerHub and start a container.  
```console
docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest 
```
5. You can now use iDEP locally by starting your web browser and enter localhost:3838 in the address bar.

Note that the Docker engine is now running in the backgroup, acting as a webserver. It works even if you restart your computer. To stop it, run these from the Terminal app: 
```console
docker stop idep 
docker rm idep
```
After stopping it, you can restart it by repeating Step 4, which also pulls the latest Docker image from DockerHub. Make sure you upgrade your iDEP image at least on a monthly basis.

Alternatively, you can also install iDEP as an R package or copy the iDEP code locally. The method is the same as above in the Windows section.

### Linux: Docker container (~10 minutes)
1. Install Docker Engine follow instructions [here](https://docs.docker.com/engine/install/). Many Linux systems have Docker installed by default. On Ubuntu, I used the following scripts.
```console
curl -fsSL https://get.docker.com -o install_docker.sh
sudo sh install_docker.sh
```
2. Pull the iDEP Docker image from DockerHub and start a container.
```console
sudo docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest 
```
3. You can now use iDEP locally by starting your web browser and enter **localhost:3838** in the address bar. If this is an server, make sure that the port 3838 is open. Then replace the localhost with the IP address.

Note that the Docker engine is now running in the backgroup, acting as a webserver. To stop it: 
```console
sudo docker stop idep 
sudo docker rm idep
```
After stopping it, you can restart it by repeating Step 2, which also pulls the latest iDEP image from DockerHub. We update it frequently, make sure you upgrade your image at least on a monthly basis.

