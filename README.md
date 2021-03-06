# verification_sdc2

## Current status
The script reads the development fits file, runs Sofia, produces a Pandas dataframe and processes the units of R.A., Dec. and size (arcsec). It also outputs a catalog that looks like the following.

```
id ra dec hi_size line_flux_integral central_freq pa i w20
79 180.2677 -29.9839 17.3036 0.1290 1108040000.0000 183.3910 26.0550 79.4879
26 179.7948 -29.9086 16.6968 0.0904 1122100000.0000 55.6169 41.1911 51.8660
24 179.8295 -30.0319 13.4996 0.0597 1122650000.0000 229.9217 35.1623 43.6766
76 180.1606 -29.9273 15.8499 0.0520 1109100000.0000 244.7785 28.7093 37.4970
93 180.1237 -30.0897 25.6036 0.0425 1106580000.0000 155.8011 22.8365 23.6625
55 180.0979 -29.8010 11.3956 0.0412 1114180000.0000 265.8551 40.3855 39.9109
107 180.0324 -30.1860 8.8770 0.0400 1100280000.0000 181.0144 42.1247 26.5400
16 180.2506 -29.9673 14.0804 0.0396 1126690000.0000 242.6396 44.2246 34.3879
....
```
#### Things to do
Basically confirm that the parameters produced are right: that the sofia definition matches what it is expected for the output catalog.


## How to use
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

### In Add_FAT branch
Rewritten the original script to only fit a single cube. Currently set for evaluation, development_small and a debug option. Need to add the large cube option.
Also now results are written in a directory of choice no longer in the main github file. And the download is only done when the cube is not already present.
Sofia is also only ran when the catalogue is not found and cubelets are produced. Added some options to sofia input and the output and input are now set dynamically.
Added a frequency to km/s convert function and converting the freq and w20 when the source finding is done in frequency.

With a trigger pyFAT can be called. This creates a fat directory in the results directory and copies the cubelets into the required per galaxy directory structure.
Also converts the cubes when in frequency.
pyFAT is available  on https://github.com/PeterKamphuis/pyFAT and can be installed through 'pip install path/to/githubclone'

It is not yet available through pip in general as it is still in development.


###To Do
Before running FAT a check needs to be done whether there is a source to be fitted.
FAT should be called internally not through os.system however currently not able to replace the command line options in call
FAT output should be compared to the initial SoFiA output and some evaluation on the accuracy of the fit.
Fix many bugs in FAT.
