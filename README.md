# iDEPGOLEM

## Description

This is the Golem shiny package of [iDEP](<https://github.com/iDEP-SDSU/idep>).
The reason for this repository is to allow for users to locally install a package version of iDEP.

## Installation

### Get the database

The data can be retrieved from this link:<http://18.235.92.206:8080/>

Once the data is unzipped, delete the data_go subfolder, which should be replaced by the most recent version from 
[here.](<https://sdsu.box.com/s/y6prfukb8pw796ya00ot2eubm4sb4svq>). This contains more demo data files.

Once all the data is retrieved then unzip the data into a base folder, then set up
an environment variable name ```IDEP_DATABASE```

How to make environment variable for each OS:

* [Windows](https://docs.oracle.com/en/database/oracle/machine-learning/oml4r/1.5.1/oread/creating-and-modifying-environment-variables-on-windows.html)
* [Mac](https://phoenixnap.com/kb/set-environment-variable-mac)
* [Linux](https://linuxize.com/post/how-to-set-and-list-environment-variables-in-linux/)

### Get the code

You need devtools(for Windows you also need RTOOLS)

```{R}
> devtools::install_github("https://github.com/espors/idepGolem")
# OR you can do the following
$ git clone https://github.com/espors/idepGolem.git
$ cd idepGolem
> devtools::install()
```

You might get the following output but it still works

```{R}
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

For some Linux OS igraph package may not install correctly look at this [issue](https://github.com/igraph/rigraph/issues/275)
and then retry the above commands.

How to run the app locally in the browser use the following command

```{R}
> idepGolem::run_app() 
Loading required package: shiny

Listening on http://127.0.0.1:7045
Browsing http://127.0.0.1:7045
```
