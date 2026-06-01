# AverageCellLoc User Guide

## Purpose

AverageCellLoc builds an average cell for _S. aureus_ cells given microscopy data. It achieves this by either averaging fluorescence or **single spot localizations**. 

Workflow: 
 - **eHooke (v1.1)** for segmentation
 - Alignement of individual cells by their major axis
 - Averaging fluorescence or spot localizations into a representative cell model. If spot models are needed **TrackMate XML files NEED to be available**. More info below

## Requirements

- eHooke 1.1 and its dependencies (Python 3.6)
- JupyterLab.
- FIJI Installation compatible with TrackMate script - [HERE](https://gist.github.com/antmsbrito/f2250a1a905457436532ee761fa6eab7)


## Quick Start

1. Prepare root folder containing your rawdata, including channel splitting.
2. Generate TrackMate XML files if spot localization will be included.
3. Create a folder to save your results. 
3. Clone or download the repo, copy the batch notebook and edit the paths and image names.
4. Run the notebook cells in order.

## Inputs

### Required

- Root folder containing one folder per experiment

```text
ROOT FOLDER
├── Experiment_1
│   ├── FoV_1
│   │   ├── Fluor_1.tif
│   │   ├── Fluor_1.xml   # optional
│   │   ├── Fluor_2.tif   # optional
│   │   ├── Fluor_2.xml   # optional
│   │   ├── Originalfile.czi  # The original image in the original file format can stay in the folder
│   │   └── Base.tif
│   ├── FoV_2
│   │   ├── Fluor_1.tif
│   │   └── Base.tif
│   └── ...
├── Experiment_2
│   └── ...
└── ...
```


### Optional

- Second fluorescence channel.
- Membrane channel.
- DNA channel for cell-cycle classification.
- TrackMate `.xml` file for spot localization.
- Pixel size value used when converting spot coordinates.



### Notes on Naming

- Keep image names consistent across all fields of view.
- Use the same base image name in every replicate folder.
- If you use a second fluorescence channel, keep its file name consistent too.
- If XML files are present, they should match the corresponding fluorescence image name.

## Notebook Setup

Before running the example notebook, fill in the following variables or cells:

- `eHooke_path`: [add the local path here]
- `root_path`: [add the root experiment folder here]
- `fluor1_name`: [add the main fluorescence file name here]
- `fluor2_name`: [add the second fluorescence file name here, if used]
- `base_name`: [add the base image file name here]
- `base_type`: [choose Phase, BF, or Membrane]
- `pxsize`: [add the pixel size here, if spot coordinates are used]
- `memb_name`: [add membrane image name here, if used]
- `dna_name`: [add DNA image name here, if used]

### Base Image Type

- `Phase`: [describe when to use this]
- `BF`: [describe when to use this]
- `Membrane`: [describe when to use this]

## Processing Pipeline

### 1. Image discovery

- The code scans the root folder for replicate subfolders.
- Each replicate folder is expected to contain the images for one field of view.
- Placeholder: describe any lab-specific folder naming rules.

### 2. Segmentation with eHooke

- The base image is loaded and masked.
- eHooke computes cell segments from the base image.
- The mask settings vary depending on the selected base image type.
- Placeholder: mention any parameter values your lab usually changes.

### 3. Cell alignment

- Each cell is aligned to its major axis.
- Fluorescence masks are rotated to match the aligned cell geometry.
- Placeholder: explain why alignment is important for averaging.

### 4. Model building

- Average-cell fluorescence models can be built for one or two channels.
- Spot-based models can be built if TrackMate XML is available.
- Cells can be filtered by spot count and by cell-cycle phase.

### 5. Optional spot localization

- TrackMate XML files are parsed to recover spot coordinates.
- Spot coordinates are converted into cell-local coordinate frames.
- Placeholder: document the exact XML export workflow your lab should use.

### 6. Optional cell-cycle classification

- If membrane and DNA images are available, cells can be classified by phase.
- The bundled model expects the classifier inputs used by the current code.
- Placeholder: add notes about microscope settings or data types that work best.

## Main Outputs

### Per-replicate outputs

- Segmentation label images.
- Aligned cell masks.
- Aligned fluorescence masks.
- Placeholder: list any files written by your notebook or downstream scripts.

### Aggregate outputs

- Average fluorescence model.
- Average spot-localization model, if XML files are provided.
- Counts of selected cells, total cells, and total spots.
- Placeholder: describe file names, plots, or saved arrays here.

## Suggested Notebook Workflow

1. Confirm the folder structure matches the expected layout.
2. Edit the notebook variables.
3. Run the import and setup cells.
4. Run the replicate loading and segmentation cells.
5. Inspect the segmentation and alignment QC outputs.
6. Build the requested model type.
7. Save or export the final result.

## Quality Control Checklist

- [ ] Are the base images being segmented correctly?
- [ ] Are the cell masks aligned in the expected orientation?
- [ ] Are fluorescence channels being loaded from the correct files?
- [ ] Do the TrackMate XML files match the image ordering?
- [ ] Are spot counts plausible for the selected population?
- [ ] Does the final average model look biologically sensible?

## Troubleshooting


## Parameter Reference

### Model selection

- `modeltype = spot`: [describe the output]
- `modeltype = average`: [describe the output]

### Cell filters

- `minspots`: [describe default and meaning]
- `maxspots`: [describe default and meaning]
- `cellcycle`: [describe valid values and meaning]

### Channel selection

- `channel = 1`: [describe channel 1]
- `channel = 2`: [describe channel 2]
