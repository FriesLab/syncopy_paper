# Syncopy Paper -- Online Materials

The demo pipeline for the [Syncopy paper](https://doi.org/10.1101/2024.04.15.589590). It uses [Syncopy](https://github.com/esi-neuroscience/syncopy) to analyze data from the [Visual Coding Project of the Allen Brain Institute](https://allensdk.readthedocs.io/en/latest/visual_coding_neuropixels.html), described in [Siegle et al., 2021](https://doi.org/10.1038/s41586-020-03171-x).


## Instructions for Reproducing the Outputs (Python and Matlab scripts)

Before you start, install git and conda if you do not have them yet, and make sure both are available as commands in your terminal. Then follow the instructions below.

### 1. Install Syncopy and Allen SDK

In the first step, we create a fresh conda environment and install both Syncopy and the Allen SDK into it:

```shell
conda create -n syncopyallen310 python=3.10
conda activate syncopyallen310
conda install -c conda-forge esi-syncopy
pip install allensdk
```

### 2. Use the script to download Allen's data

Now, you are ready to download and run the scripts contained in this repository. First, you will have to run the script that downloads the Allen data. The script `read_allen.py` does this for you. It will take some time to download all the data of one sample session (> 10 GB). Once downloaded, future invocations of the script will use the local data instead of re-downloading it. When all data has been downloaded, the script selects the parts of the data that are used in the analyses pipeline and saves these parts in both Fieldtrip and Syncopy formats.


```shell
git clone https://github.com/frieslab/syncopy_paper
cd syncopy_paper/code/Python/
python ./read_allen.py
```

Be patient while the data is downloaded, selected, and saved in the required formats. All the later scripts require the files created by this script, so check that it was completed without errors. If everything went well, you should have the following files (in the given sub directories) in your working directory, i.e. in `syncopy_paper/code/Python/`:


* file `allen_lfp.spy/allen_lfp.analog`, size ~ 153 MB, md5sum a0611868b8d02c786a8e43b0c9735e98
* file `allen_spike.spy/allen_spike.spike`, size ~ 2 MB, md5sum 90c6e13e6daec2b5e638ac1c74dc0188
* file `allen_FT.mat`, size ~ 215 MB, md5sum 15a6f300b2ecf6caefbf869395bab74c


### 3. Run the pipeline you want

You can now run any of the Python pipeline scripts in the `code/Python/` directory. The following three scripts are available:

* Raster_PSTH.py
* PowerSpectrum.py
* ConnectivityAnalysis.py

E.g., to run the first one, type `python ./Raster_PSTH.py`.

This will run the analysis pipeline and create the output plots in PDF and PNG format in the current working directory.


### 4. Matlab code
To run the Matlab script, first you need to download the [Fieldtrip toolbox](https://www.fieldtriptoolbox.org/).
Next, go to `code/Matlab` in Matlab, open the Matlab script ```Matlab_power_connectivity.m``` and set the path to the Fieldtrip toolbox you downloaded in the Matlab script. 
Then, you can run the script in Matlab.
