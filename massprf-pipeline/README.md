# massprf-pipeline
pipeline for genome wide MASSPRF analysis

## Installation
Note: These instructions are by no means exhaustive and will change based on your individual shell configurations and environment.  In general, this instruction (and the pipeline itself) assumes a unix environment running bash.  Advanced users w/ different configurations should be able to install without these instructions.
#### 0) Install massprf and massprf_preprocess.  You can clone this repository and read within for build instructions. https://github.com/zimingz/MASSPRF_10July2016

It is important to note that you may need to build both massprf and massprf_preprocess independently.  The final executable for massprf will be in ./bin, while the massprf_preprocess executable will be in ./MASSPRF_preprocessing_08July2016.

##### 0.1) Build symlinks to massprf and massprf_preprocess

Get the absolute path of your compiled massprf & MASSPRF_preprocess.  These will be something like:

~/PATH/TO/MASSPRF/MASSPRF_10July2016/bin/massprf
~/PATH/TO/MASSPRF/MASSPRF_10July2016/MASSPRF_preprocessing_08July2016/MASSPRF_preprocess

Create a custom symlinks folder in your home directory if you haven't already.
`mkdir ~/symbolics`

change to that directory:

`cd ~/symbolics`

Create links to massprf & massprf_preprocess in that directory

`ln -s ~/PATH/TO/MASSPRF/MASSPRF_10July2016/bin/massprf massprf`

`ln -s ~/PATH/TO/MASSPRF/MASSPRF_10July2016/MASSPRF_preprocessing_08July2016/MASSPRF_preprocess massprf_preprocess`

Add the symbolic link folder to your $PATH in ~/.bash_profile

`vim ~/.bash_profile`

Scroll to the line that says something like 

`PATH=$PATH:OTHERPATHS`

Append a colon and add the path to your symbolic links, such that it looks like this:

`PATH=$PATH:OTHERPATHS:$HOME/symbolics`


#### 1) Install conda package manager:

`pip install conda`

If you are working on a cluster, you may need to install miniconda directly from the package due to file permissioning/pip being unavailable.  In this case, download the file from https://conda.io/miniconda.html , run  it on the cluster, and make sure that the installation directory is added to your path via ~/.bashrc.  After installing miniconda, use `source ~/.bashrc` to source the path, then type `python` and check that your python version is >3.5 and is produced by 'Continuum Analytics'

#### 2) Make sure your conda is up to date by running:

`conda update conda`

#### 3) Update python to 3.5, or (optionally) create a python 3.5 virtual environment:

`conda update python`

#### 4) Get package dependencies

Add Bioconda channel to Conda: 

`conda config --add channels bioconda`

Pandas:

`conda install pandas`

Biopython:

`conda install biopython`

Gffutils: 

`conda install gffutils`

pyvcf: 

`conda install pyvcf`

#### 5) Clone this repository

`git clone https://github.com/Townsend-Lab-Yale/massprf-pipeline.git`


# Usage:

To run example, see jobs.list
