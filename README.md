# AverageCellLoc
## Averaging fluorescence in _S. aureus_ cells

This is is a tool to build an heatmap of the fluorescence of an average _S. aureus_ cells. The code heavily relies on [eHooke's][1] (version 1.1) and optionally [StarDist][2].

Given a set of microscopy images the tool automatically runs eHooke in order to segment cells, orients each cell according to their
major axis and averages the fluorescence in order to build a model cell. It requires per FoV a base image that eHooke can segment (phase, brightfield or membrane dye) and it works with 1 or 2 fluorescence channels.

If spot detection was done using [TrackMate][3] and the output .xml file is provided it can
also build an heatmap of the average spot localization. Spot detection HAS to be done using the following gist - https://gist.github.com/antmsbrito/f2250a1a905457436532ee761fa6eab7 (otherwise the order of the XML nodes dont match)


## Requirements
####  Built in Python 3.6
#### - eHooke 1.1 and its dependencies 

## Usage instructions

#### - Organize image data. 
```
ROOT FOLDER
│
├─── Experiment #1 
│    │
│    ├───FoV_1
│    │   ├──   Fluor_1.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   ├──   Fluor_2.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   └──   Base.tif
│    ├───FoV_1
│    │   ├──   Fluor_1.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   ├──   Fluor_2.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   └──   Base.tif
│    └─── ...
├─── Experiment #2
│    │
│    ├───FoV_1
│    │   ├──   Fluor_1.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   ├──   Fluor_2.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   └──   Base.tif
│    ├───FoV_1
│    │   ├──   Fluor_1.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   ├──   Fluor_2.tif
│    │   ├──   Fluor_1.xml (optional)
│    │   └──   Base.tif
│    └─── ...
...

```
#### - Run the jupyter notebook 
You have to provide the following information:
 1. The path to your local eHooke folder
 2. The path to the root folder where the images are (see above for example)
 3. The file name of the fluorescence images (without .tif)
 4. The file name of the base images (without .tif)
 5. The type of base image provided

You have to edit the example jupyter notebook file according to your needs.


[1]: https://github.com/BacterialCellBiologyLab/eHooke/releases/tag/v1.1.0 
[2]: https://doi.org/10.1007/978-3-030-00934-2_30
[3]: https://doi.org/10.1038/s41592-022-01507-1

