# Sites of Metabolism (SoMs) Data Analysis
Code used in preprint [Metabolite Identification Data in Drug Discovery: Site-of-Metabolism Annotation, Analysis, and Exploration for Machine Learning](https://doi.org/10.26434/chemrxiv-2025-6mtq8) 

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



