# CRAMPON_v1.0: pre-post processing module and plotting scripts (GMD manuscript version)
---

## General description
---

CRAMPON (CRocus with AssiMilation of snowPack ObservatioNs) is an ensemble data assimilation system assimilating snowpack observations with a Particle Filter.

The present package "crampon" is a python toolbox to launch, pre/post-process and plot CRAMPON simulations on Météo-France HPC system, in link with the publication manuscript submitted to _Geoscientific Model Development_ (GMD):
>CRAMPON : a Particle Filter to assimilate snowpack observations
in a spatialised framework,
Bertrand Cluzet, Matthieu Lafaysse, Emmanuel Cosme, Clement Albergel,
Louis-François Meunier, and Marie Dumont, (submitted) **(doi://...)**

## Main Features
---
CRAMPON can be run without the present package **crampon**, but only on Météo-France HPC environment, using **snowtools** and **vortex** packages (see **snowtools** [CRAMPON user doc](https://opensource.umr-cnrm.fr/projects/snowtools_git/wiki/CRAMPON_user_doc)). **crampon** completes this packages by providing additional features to pre-post-process and launch CRAMPON experiments _in any environment_ :
- launching (numbers of) openloop and synthetic data assimilation experiments.
- generating synthetic observations from openloop experiments, masking them (```CramponObs```).
- post-processing thousands of netCDF outut files from SURFEX-Crocus into a handful of easily-readable pickle files (pickle2 and 3 handled).
- loading these pickle files into python objects (```CramponPp```) making easy to manipulate/plot.
- several plotting methods tailored to the model's semi-distributed geometry (```Pie``` in particular).
- CRAMPON's Particle Filter can be launched standalone locally for development/tests purposes (```CramponPf```).
- This package also includes tools to locally launch parallelized CRAMPON sequences (```CramponParallel```) _(beta version)_. This can be useful for testing and developping on small domains/ensemble, but deployment on an HPC environment is highly recommended for bigger application.

### In addition:
- Notebooks producing the manuscript figures are provided for reproducibility (The companion dataset can be downloaded at...) and to serve as examples.
- additional examples on how to use this package can be found in the ```examples``` repertory.


## Dependencies
---
- **Python libraries**: The code is developped on latest stable version of python3, and has been tested with python 3.5.2. [Required libraries](link_to_requirements.txt): matplotlib, numpy, pandas, scipy, netCDF4, seaborn, pytest, configobj. Installation inside a conda virtual environment is recommended.
- **snowtools**: open-source python package used to pre/post process SURFEX-Crocus outputs and launch simulations on Météo-France HPC environment. Used for operational snowpack modelling at Météo-France. The present package is highly inspired on snowtools and will be included in it in the future. Please carefully check snowtools documentation (in particular the [wiki](https://opensource.umr-cnrm.fr/projects/snowtools_git/wiki/Wiki) and its [CRAMPON user doc](https://opensource.umr-cnrm.fr/projects/snowtools_git/wiki/CRAMPON_user_doc)) as it includes detailed information on the installation procedure and the simulation environment.
Download, documentation and instructions : https://opensource.umr-cnrm.fr/projects/snowtools_git/ Manuscript version: tag *CRAMPON_v1.0*
- **SURFEX** : open-source Land Surface Model, including FORTRAN sources of CRAMPON snowpack model (SURFEX-ISBA-Crocus) and the Particle Filter.
Download and instructions : https://opensource.umr-cnrm.fr/projects/surfex_git2/ Manuscript version: tag *CRAMPON_v1.0*

- **vortex** (optional): python package gathering all environment-specific codes of Météo-France modelling systems relative to its HPC computing system. Can not be applied in other HPC environment. Install is not mandatory. A tool to locally launch CRAMPON in parallel is provided in the present package as an alternative for low computational cost applications _(beta mode)_.
Download and instructions : https://opensource.umr-cnrm.fr/projects/vortex/ Manuscript version : tag *CRAMPON_v1.0*



## Install
---

Once you installed the required dependencies (see above) you can clone this repository.

In your .bashrc or equivalent, add **crampon** root to your $PYTHONPATH:
```bash
export PYTHONPATH=<your_path_to_crampon>:$PYTHONPATH
```

As this package closely follows vortex's file structure, it is necessary to put all your archive/simulation in a specific ```path```. In your .bashrc or equivalent set the following environment variable:
```bash 
export CRAMPONPATH=<path>
```



## Getting started
---

First of all, thouroughly read the documentation of CRAMPON within the snowtools module :[CRAMPON user doc](https://opensource.umr-cnrm.fr/projects/snowtools_git/wiki/CRAMPON_user_doc), in order to get familiar with the software environment. 
- Basically, the only inputs you needs are a meteorological forcing and a set of observations. Forcings from SAFRAN reanalyses (_Vernay et al., (in prep)_) over French mountain ranges are available at https://doi.org/10.25326/37. Observations must be converted to the daily format specified in the doc. 
- Use this forcing to generate a *spinup* (PGD.nc and PREP.nc files) using **snowtools** _s2m command_ (example in pproc_scripts/spinup.py) and an appropriate namelist (basic examples in snowtools_git/DATA/).
- Apply stochastic perturbations on the forcing to generate an *ensemble of forcings* (snowtools/tools/makeForcingEnsemble.py or snowtools_git/tools/job_gener_pert_forcings.py)
- Prepare the assimilation sequence : create a configuration file with the assimilation dates and the ESCROC member ids. choose your PF configuration (to be directly written in the namelist).
- launch the assimilation sequence (or an openloop) on Meteo-France HPC system (inspiring on pproc_scripts/multilaunch.py) or locally (inspiring on examples/launch_parallel_local.py).
- the [experiment/archive map](link_to_thexmind_map) tries to give you an overview of the experiment and archive file structure. 
- post process the archive into pickle files (using CramponPp.py, inspiring on examples/launch_parallel_local.py or pproc_scripts/multipp.py)
- once an experiment has succeeded, you can also play around with the Particle Filter, applying it on the background PREPS of specific dates (using CramponPf.py, inspiring on examples/local_pf_run.py).
- enjoy :)

## License information
---

See the [license](https://github.com/smrt-model/smrt/blob/master/LICENSE) (copy from smrt, replace with real license file). file for terms & conditions for usage, and a DISCLAIMER OF ALL WARRANTIES.

DISCLAIMER: This version of CRAMPON is under peer review. Please use this software with caution, ask for assistance if needed, and let us know any feedback you may have.

Copyright (c) 2020 Bertrand Cluzet, Matthieu Lafaysse, Marie Dumont.

## Other Contributions
---

- César Deschamps-Berger (plotting routines, not used in the manuscript).


