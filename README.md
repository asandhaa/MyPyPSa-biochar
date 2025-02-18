# MyPyPSa-biochar

## Contents

- [Introduction](#introduction)
- [Repository structure](#repository-structure)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
  
## Introduction
MyPyPSA-biochar, a myopic optimization model developed to represent the German energy system with a detailed mapping of the electricity sector, on a highly disaggregated level, spatially and temporally, with regional differences and investment limitations. Pyrolysis is implemented as a power generator and a negative emission technology.

MyPyPSA-Ger was developed by [EEW group](https://ines.hs-offenburg.de/forschung/energiesysteme-und-energiewirtschaft) at [Hochschule Offenburg](https://www.hs-offenburg.de/). 
MyPyPSA-biochar model is built using the Modeling Framework [MyPyPSA-Ger](https://github.com/AnasAbuzayed/MyPyPSA-Ger/), and upon [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur).

The model is described in the paper [Pyrolysis as a Strategic Element in Energy System Transformation to Achieve Net Zero Emissions](Include Link) 




## Repository structure

The model requires various input data files and network configurations to run simulations effectively. Below is an overview of the directory structure and the contents relevant to the energy system model.


### Code
In the main directory, you will find several Python files that contain the model's code. Among these:
- **config.yaml**: This is where basic settings are defined and the initial conditions of the network before running simulations.
- **Model.py**: This is the coordinator file that orchestrates the execution of the entire energy model by calling functions from other files.
- **Other Python Files**: These contain specific functions and modules that support the overall functionality of the model.

### Data
In the main directory, you will find several folders model's data. Among these:

- The `data/` folder contains essential datasets that the model uses for simulations:
  - Installed power plant capacities until 2019
  - costs
  - co2-limits
- The `1_Regionality/`folder contains the regional biomass potential for pyrolysis
- The `scripts/`folder contains additional codes for solving
- The `starting networks/` folder holds pre-configured network files that define the initial structure and components of the energy system model. These files allow you to start simulations of the BAU (network without pyrolysis) or Pyrolysis scenario (network with pyrolysis)
- The `elec_s_37_ecb_lcopt_Co2L-1H-Ep/` folder serves as container for the networks for all years and already contains 37 cluster related data generated before the myopic optimisation

### 

## Installation 

### Clone the Repository 

/some/other/path % cd /some/path/without/spaces

/some/path/without/spaces % git clone https://github.com/asandhaa/MyPyPSA-biochar.git

### Install the Library

% cd MyPyPSA-biochar

% conda create --name MyPyPSA-biochar --file environment.yml

### Download a Solver of your choice

## Usage
### 1. BAU and Pyrolysis scenario
To run either the BAU or the Pyrolysis scenario, only the correct starting network from the folder `starting networks/` needs to be copied into the main folder and the prefix [BAU scenario] or [Pyrolysis scenario] needs to be deleted.

### 2. Scenario settings
Scenario setting can be adapted:

#### General

- CO2 price in `data/co2_price.csv`
- CO2 limit in `data/co2limits.csv`
- Yearly increase in total electricity demand in `Model.py` in function Myopic.update_load(n,1.01). 1.01 is +1%/year, 1.02 would be +2%/year, ...
- Yearly installation limit in `data/agg_p_nom_minmax.csv`

#### Pyrolysis related parameter
- CAPEX of pyrolysis in `data/capex_pyro.xlsx`
- OPEX directly in `Model.py`
- National biomass potential for pyroylysis in `config.yaml`
- Yearly installation limit for pyrolysis in `pyrolysis.py` in function add_CCL_constraints_pyro()

### 3. Run MyPyPSA-biochar model

```bash
conda activate MyPyPSA-biochar
python Model.py
```
% enter data network name
```bash
elec_s_37_ecb_lcopt_Co2L-1H-Ep-CCL.nc
```
% enter regional potential value
```bash
2500
```

The results of this model will be saved in the folder `/elec_s_37_ecb_lcopt_Co2L-1H-Ep` within the main repository.
## License
This work is licensed under CC BY 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by/4.0/
