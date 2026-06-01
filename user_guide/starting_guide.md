# AverageCellLoc User Guide

## Purpose

AverageCellLoc builds an average cell for _S. aureus_ cells given microscopy data. It achieves this by either averaging fluorescence or **single spot localizations**. 

Workflow: 
 - **eHooke (v1.1)** for segmentation
 - Alignement of individual cells by their major axis
 - Averaging fluorescence or spot localizations into a representative cell model. If spot models are needed **TrackMate XML files NEED to be available**. More info below

## Inputs

### Required

- Root folder containing one folder per experiment
- Folder to save results to. 
- Download eHooke 1.1, AverageCellLoc and create a suitable python 3.6 environment
- FIJI and Trackmate (optional)
- The following scripts for preprocessing
    - Splitting and reorganizing raw data into the expected format - [HERE](https://gist.github.com/antmsbrito/7d9602ad6da8c84a1a37b9765bce5f58)
    - Exporting TrackMate XML files - [HERE](https://gist.github.com/antmsbrito/f2250a1a905457436532ee761fa6eab7)

## STEP #1 Set up your root folder

Root folder should contain one folder per experiment. Each experiment folder should contain the raw data for all FoVs of that experiment. The expected folder structure is as follows:

```root_folder/
├── experiment_1/
│   ├── FoV1.czi
│   ├── FoV2.czi
│   └── ...
├── experiment_2/
│   ├── FoV1.czi
│   ├── FoV2.czi
│   └── ...
└── ...
``` 

Open FIJI and use the following [script](https://gist.github.com/antmsbrito/7d9602ad6da8c84a1a37b9765bce5f58) to split and reorganize the raw data into the expected format. This script assumes that your raw data is in CZI format. The script will create separate folders for each experiment and save the split images in the appropriate format. 

IMPORTANT: CHANGE LINE 45 AND 49 OF THE SCRIPT TO MATCH YOUR IMAGE AND FOLDER STRUCTURE. 

For example if your raw data is in a folder called `rawdata` and your images have 2 channels named `channel1` and `channel2`, you would change line 45 to:

```python
rootpath = r"rawdata" 
```
And line 49 to:

```python
		name = ['channel1','channel2']
``` 

To use it, just drag and drop it into FIJI and click run. The folder structure now should look like this:

```root_folder/
├── experiment_1/
│   ├── 1
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   └── FoV1.czi
│   ├── 2
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   └── FoV2.czi
│   └── ...
├── experiment_2/
│   ├── 1
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   └── FoV1.czi
│   ├── 2
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   └── FoV2.czi
│   └── ...
└── ...
``` 


## STEP #2 Export TrackMate XML files (optional)

This step is optional but required if you want to build spot-based models. If you only want to build fluorescence-based models, you can skip it.

Open FIJI and use the following [script](https://gist.github.com/antmsbrito/bc824dc9456223e27ffb13e7a1bc0bd2) to export TrackMate XML files. This script assumes that your images are in the format described in STEP #1. The script will create a TrackMate XML file for each FoV **NOT NAMED PHASE OR LABELS** and save it in the appropriate folder.

IMPORTANT: CHANGE LINE 35, 48-53 AND 93 TO MATCH THE DESIRED SETTINGS. 

Line 35 should be changed to match the pixel size in microns of your images: 

```python
    IJ.run(imp, "Properties...", "channels=1 slices=1 frames=1 pixel_width=0.0645 pixel_height=0.0645 voxel_depth=1.0000");
```

Line 47-53 should be changed to match the settings you want to use for spot detection. For example, if you want to use a LoG detector with a **radius** of 0.2 microns you would change it to:

```python
settings.detectorSettings = {
        'DO_SUBPIXEL_LOCALIZATION': True,
        'RADIUS': 0.2,  # CHANGE IF NEEDED
        'TARGET_CHANNEL': 1,
        'THRESHOLD': 1.,  
        'DO_MEDIAN_FILTERING': False,
    }
``` 

Line 93 should be changed to your rootfolder path. 


To use it, just drag and drop it into FIJI and click run. The folder structure now should look like this:

```root_folder/
├── experiment_1/
│   ├── 1
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   ├── FoV1.czi
│   │   ├── channel1.xml
│   │   └── channel2.xml
│   ├── 2
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   ├── FoV2.czi
│   │   ├── channel1.xml
│   │   └── channel2.xml
│   └── ...
├── experiment_2/
│   ├── 1
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   ├── FoV1.czi
│   │   ├── channel1.xml
│   │   └── channel2.xml
│   ├── 2
│   │   ├── channel1.tif
│   │   ├── channel2.tif
│   │   ├── FoV2.czi
│   │   ├── channel1.xml
│   │   └── channel2.xml
│   └── ...
└── ...
```



## STEP 3 - Run the notebook

Make sure you have the required dependencies installed. 

eHooke 1.1 is required and should be working. Confirm that you also have Jupyter Lab installed and working. In the terminal, you can check if Jupyter Lab is installed by running:

````
pip install jupyterlab
````

Download AverageCellLoc and make a copy of the `RunCA_Batch_example.ipynb` notebook and rename it to something you prefer. 

Open the notebook and run it. You can run the notebook in Jupyter Lab by navigating to the folder where the notebook is located and running:

````bash
jupyter lab RunCA_Batch_example.ipynb
````

Before running the notebook, read the instructions and fill in the following variables or cells:

- `ehooke_path`: [add the path to eHooke 1.1 here]
- `root_folder`: [add the path to the root folder here]
- `result_folder`: [add the path to results folder here]

- `base`: [name of the filename to be used as the base image (img to be segmented)]
- `membrane`: [name of the filename containing membrane signal, if available - ONLY USED FOR CELL-CYCLE CLASSIFICATION]
- `basetype`: [type of the base image, either 'Phase', 'BF' or 'Membrane' - important for segmentation]
- `dnaname`: [name of the filename containing DNA signal, if available - ONLY USED FOR CELL-CYCLE CLASSIFICATION]
- `pxsize`: [pixel size in microns of the images] 

- `fluor1`: [name of the filename containing the first fluorescence signal]
- `fluor2`: [name of the filename containing the second fluorescence signal, if available, otherwise set to 'none']


### The following is also VERY IMPORTANT to check: 

- `model.py` : Check eHooke settings and make sure they are suitable for your data. Although most are default some are important to check such as channel alignment and mask dilation. 
