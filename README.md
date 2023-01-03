# Cell Averager
## Averaging fluorescence in _S. aureus_ cells

Cell Averager is a tool to build an heatmap of the fluorescence of an average _S. aureus_ cells. The code heavily relies on [eHooke's][1] capabilities to segment cells.

Given a set of microscopy images the tool automatically runs eHooke in order to segment cells, orients each cell according to their
major axis and averages the fluorescence in order to build a model cell. It requires per FoV a base image that eHooke can segment (phase, brightfield or membrane dye) and it works with 1 or 2 fluorescence channels.

If spot detection was done using [TrackMate][2] and the output .xml file is provided it can
also build an heatmap of the average spot localization.

This tool (without spot detection support) is implemented in this [eHooke branch][1].

## Requirements
####  Built in Python 3.6
#### - eHooke and its dependencies (https://github.com/antmsbrito/eHooke_1.0)

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
#### - Run `main.py` from eHooke's environment and follow the instructions in the GUI
You have to provide the following information:
 1. The path to your local eHooke folder
 2. The path to the root folder where the images are (see above for example)
 3. The file name of the fluorescence images (without .tif)
 4. The file name of the base images (without .tif)
 5. The type of base image provided
 6. Press ENTER!

## Acknowledgements

This tool relies heavily on a fork of eHooke. eHooke is a tool originally developed by [Bruno Saraiva][3]. The original
eHooke publication can be found [here][4]. 

Furthermore, the current implementation of this tool is very heavily inspired by an earlier tool developed by [Bruno][5].
Specifically the `CellAligner` and `CellModeler` classes are refactors of classes he wrote. Therefore, the original idea
behind this tool can be attributed to him.

This tool has a repository of his own in my personal page with his permission. 

[1]: https://github.com/antmsbrito/eHooke_1.0
[2]: doi:10.1038/s41592-022-01507-1
[3]: https://github.com/brunomsaraiva
[4]: doi:10.1017/S2633903X21000027
[5]: https://github.com/brunomsaraiva/CellAverager