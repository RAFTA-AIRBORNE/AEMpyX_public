# Tools for processing, inversion, and interpretation of Airborne ElectroMagnetics (AEM)

This is the public repository for the  Extended Airborne Electromagnetic Inversion toolbox "AEMpyX". **It will be successively populated in the near future.** 

The "AEMpyX" toolbox was originally developed at DIAS starting with the project "Spatially constrained Bayesian inversion of frequency- and time-domain electromagnetic data from the Tellus projects" (2015-sc-004), followed by "RAFTA: Resolution Analyses for Frequency- and Time-Domain Airborne Electromagnetic Data of the Irish Tellus Programme" (2020-sc-049), both funded by the Geological Survey of Ireland GSI. It is distributed under the GNU GENERAL PUBLIC LICENSE Version 3. Currently it is expanded to include more algorithms (in particular for uncertainty quantification, UQ), platforms, and a branch for making it useful for magnetotelluric (MT) modelling and inversion. 

AEMpyX currently works fully under linux operating systems, but mostly also under  windows (here short for Windows 10). There are of course changes in the installation procedures, as mentioned below. 

Under linux, get your working copy directly via _git_ from the command line. In windows, _git_ functionality is available, once Ana/Miniconda is installed (see below). There, open the included terminal, and clone the repository with the same line as in linux: 

_git clone https://github.com/RAFTA-AIRBORNE/AEMpyX_public.git_

The created local repository AEMpyX contains the following subdirectories:

- 	**environment**
	Contains conda environment description files, and some useful helper files for working within the 
	conda environment. For linux and windows you can find special versions in the repective directories. 
	The current AEMpyX environment contains some packages which are not strictly necessary for running 
	aempy, but useful for related geoscientific work. 
 	
- 	**aempy/scripts**
 	Contains a collection of scripts for processing, visualization, and other tasks releted to one-dimensional inversion of 
 	AEM data;. Also included is a workflow for static shift correction of MT observations (work in progress).   

- 	**aempy/tutorial**
 	Contains python scripts and jupyter notebooks for the most important scripts, explaining typical AEM preprocessing, 
	visualization, one-dimensional inversion, and uncertainty quantification workflows, using the toolbox. 

-	**aempy/modules**
 	Contains the modules _aesys.py_, _prep.py_, _post.py_, _viz.py_, _inverse.py_, _alg.py_, and _util.py_, which 
	are called from the Python scripts in **aempy/scripts** and **aempy/tutorial**, accomplishing 
	different tasks of AEM inversion. It also contains the _core1d.so_ (or _core1d.pyd_ whem using windows) module, 
	once compiled from the sources in **aempy/core1d**.
 	
- 	**aempy/core1d**
	This directory contains the Fortran 90 source code for the computational core run by the Python toolbox. 
	The numerics is derived from the AMIRA/CSIRO AirBeo software. Currently it contains working wrappers for 
	the two systems used in Tellus: AEM05, and GENESIS, with AEM95 (used in the earlier Northern Ireland 
	surveys) available, and  TEMPEST and GEOTEM under development. This directory also contains makefiles 
	for both operating systems, named _Makefile_linux_ or _Makefile_windows_, respectively.   

-	**aempy/util**
 	Contains  useful helper shell scripts etc. 
	
-	**aempy/publish**
	This is the place where the latest report, and related publications, can be found.

-	**aempy/info**
 	Documentation for python (including the most important extensions, numpy,
	scipy, and matplotlib), and other useful tools.

This version will run under Python 3.9+. However, moving to 3.10/11 is not yet encouraged because they are still missing some important (though not essential) packages. To install the python environment in any Linux environment (e.g. Ubuntu, SuSE), you need to do the following:


(1) Download the latest Anaconda or Miniconda version (https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html), and install by running the downloaded bash script with:  

_bash Miniconda3-latest-Linux-x86_64.sh_

For windows just execute the downloaded executable, _Miniconda3-latest-Windows-x86_64.exe_. In order to make updates secure and avoid inconsistencies, copy _.condarc_ from AEMpyX/environment to your home directory. As the Miniconda installer is not updated very frequently, it is useful  to run the following within Anaconda:

_conda update conda_

_conda update --all_


Do this regularly to keep everything consistent! 


(2) Create an appropriate conda environment (including the necessary prerequisites) from the files AEMpyX.yml or AEMpyX.txt found in the AEMpyX/environment  directories by entering

_conda env create -f AEMpyX.yml_

or

_conda create --name AEMpyX --file AEMpyX.txt_

in the command window (powershell under windows).

This will set up a Python 3.9 (see https://docs.python.org/3.9/) environment with all dependencies for aempy. Don't forget to update also the used environment regularly, using _conda update --name AEMpyX --all_! 

There is a replacement for _conda_, called _mamba_ (see https://github.com/mamba-org/mamba), which is not only considerably faster, but also better in keeping the environments consistent. It can be installed via _conda_ (i. e., _conda install mamba -c conda-forge_), and has practically the same syntax as the original package manager.  


(3) Activate this environment by:

_conda activate AEMpyX_


(4) Within this envionment you now need to compile the aempy core, which is written in Fortran 90, and thus needs. For this purpose, _f2py_ (part of _numpy_), and the required  compilers have been included in the EM environments. To compile and install the numerical core, go to the _core1d_ directory, and enter 

_make -f Makefile_linux_  or  _make -f Makefile_windows_ 

respectively. If this is succesful, a dynamical library, _core1d.so_ (_core1d.pyd_ under windows), should be in the _modules_ subdirectory. 


(5) In order to reproduce identical behavior of matplotlib, you should copy the included  _matplotlibrc_ file to the appropriate directory. Under Linux (Ubuntu), this should be : _$HOME/.config/matplotlib/matplotlibrc_. Pertinent changes should be made there, or have to be made within the scripts/modules using the _mpl.rcParams[name]=value_ mechanism. 


(6) For running aempy scripts, we have defined two environmental variable, _AEMPYX_ROOT_ and _AEMPYX_DATA_. These point to the place where AEMpyX is installed, and where you keep your AEM data, respectively. Keeping to this scheme makes life much easier if more than one person work on the tools.

In linux you can set them in your _.bashrc_ file. Example: 

_export AEMPYX_ROOT='${HOME}/AEMpyX/'_	

_export AEMPYX_DATA='${HOME}/AEM_Data/'_

Under windows, you should use the system settings dialogue to do so. 


(7) Finally, the remaining open source toolboxes you want to use need to be installed, either via the _conda_ or _mamba_ framework, or the _pip_ command. Using _pip_ within a _conda/mamba_ environment should be kept to a minimum, as it can lead to serious inconsistencies. You may need to re-install via _pip_ again after an update of the _conda_ environment. 


(8) Once in the activated conda environment _AEMpyX_, there are several ways to start python scripts or jupyter notebooks (_https://jupyter.org/_). 

For getting started with python **scripts**, we suggest to use the _spyder_ IDE (_https://www.spyder-ide.org_/), which is already installed within the _AEMpyX_ environment. It has been developed for easy development of python software, including visualisation with _matplotlib_ and derived packages. Current versions (5.X) also allow to work with other languages as _JULIA_ or _R_, or with _jupyter notebooks_ by installing the appropriate plugin. 

However, as we have defined the environmental variables in step (6), the scripts can be run from anywhere in the system, either from the activated _AEMpyX_ environment, or using _conda run_ (see, _https://docs.conda.io/projects/conda/en/latest/commands/run.html_) from an initialized conda (i.e., base environment). Running python scripts directly from the command line can be done in different ways:

_python3 -u mypythonscript.py > output.log_, 

if you are in AEMpyX. If the magic first line _#!/usr/bin/env python3_ exists in the skript, it can be called as 

_mypythonscript.py > output.log_.

If anywhere in an conda environment, use 

_conda run_ -n AEMpyX mypythonscript.py

The usual way to work with **notebooks** is with your favorite web browser, or specialized software as _jupyterlab_. If using your browser (as set in your system as default) you can simply use the classical interface

_jupyter notebook mynotebook.ipynb_ 

or the new one

_jupyter lab mynotebook.ipynb_

Both calls will open a new browser window, in which you can edit and run the notebook. The _jupiterlab_ is the future interface for notebooks, and there are many options not available with the classical call (_https://jupyterlab.readthedocs.io_). Remote clusters often offer a specialized server _JupyterHub_) to develop and run notebooks. 

Enjoy, read the docs, but please keep in mind that this is an experimental software, and may contain errors. Use at your own risk! However, we will frequently update the repository correcting bugs, and adding additional functionality.   

D. Kiyan & V. Rath

May 5, 2023
