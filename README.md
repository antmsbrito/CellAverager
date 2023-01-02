# Cell Averager
## Averaging fluorescence in _S. aureus_ cells

Cell Averager is a tool to build an heatmap of the fluorescence of an average _S. aureus_ cells. The code heavily relies on [eHooke's][1] capabilities to segment cells.

Given a set of microscopy images the tool automatically runs eHooke in order to segment cells, orients each cell according to their
major axis and averages the fluorescence in order to build a model cell. It requires per FoV a base image that eHooke can segment (phase, brightfield or membrane dye) and it works with 1 or 2 fluorescence channels.

If spot detection was done using [TrackMate][2] and the output .xml file is provided it can
also build an heatmap of the average spot localization.

This tool (without spot detection support) is implemented in [eHooke][1].


## Requirements
####  Built in Python 3.6
#### - eHooke (https://github.com/antmsbrito/eHooke_1.0)

## Usage instructions

#### - Organize image data. 
```
Experiment 
│
└───FoV_1
│   ├──   Fluor_1.tif
│   ├──   Fluor_1.xml (optional)
│   ├──   Fluor_2.tif
│   ├──   Fluor_1.xml (optional)
│   └──   Base.tif
└───FoV_1
│   ├──   Fluor_1.tif
│   ├──   Fluor_1.xml (optional)
│   ├──   Fluor_2.tif
│   ├──   Fluor_1.xml (optional)
│   └──   Base.tif
└─── ...

```
#### - 

[1]: https://github.com/antmsbrito/eHooke_1.0
[2]: doi:10.1038/s41592-022-01507-1