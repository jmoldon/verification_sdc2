# verification_sdc2

Here are the commands to run the analysis:
```
git clone https://github.com/jmoldon/verification_sdc2.git
cd verification_sdc2/
conda env create -f environment.yml
conda activate sdc2-ver
python scripts/analysis.py
```
If you don't have conda installed, please follow these commands:

For Linux
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p conda-install
source conda-install/etc/profile.d/conda.sh
```

For MacOS, Sofia-2 does not seem to be available as a conda package. You will need to follow the instructions here: https://github.com/SoFiA-Admin/SoFiA-2. So, if you are using Mac, you will need to remove the sofia-2 dependence from the `environment.yml` file
