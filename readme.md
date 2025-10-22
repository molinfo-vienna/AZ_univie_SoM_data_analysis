# Sites of Metabolism (SoMs) Data Analysis

[![DOI](https://img.shields.io/badge/DOI-10.1021/acs.molpharmaceut.5c00740-blue)](https://doi.org/10.1021/acs.molpharmaceut.5c00740)  

Code used in [Metabolite Identification Data in Drug Discovery, Part 2: Site-of-Metabolism Annotation, Analysis, and Exploration for Machine Learning](https://doi.org/10.1021/acs.molpharmaceut.5c00740)

Examplified by [AstraZeneca site-of-metabolism (SoM) data for 120 compounds](https://zenodo.org/records/15458631) and approved drugs from [DrugBank](https://go.drugbank.com/). 

## Environment Setup


**Create the environment from `environment.yml`**
   ```bash
   conda env create -f environment.yml
   conda activate az_som_py310
   ```


**Alternative: Create a virtual environment and install dependencies**

   ```bash
   conda create -n som_env_py310 python=3.10
   conda activate som_env_py310
   conda install -c conda-forge rdkit
   pip install -r requirements.txt
   ```


## Notes
- Data files are expected in the `data/` directory.
- Outputs are written to the `output/` directory.
- Active learning results generated using code from [FAME.AL](https://github.com/molinfo-vienna/FAME.AL)



