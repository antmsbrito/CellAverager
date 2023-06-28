
import os
import sys
import tkinter as tk

from GUI import MainGUI

import numpy as np
import matplotlib
from skimage.transform import resize
from skimage.external import tifffile
from skimage.io import imsave
from skimage.color import gray2rgb

from cellmodeler import CellModeler
from cellaligner import CellAligner


def run_CellAverager(base, fluor, basetype, membrane):
    """
    Runs ehooke segmentation on two images and saves the label image
    Aligns cells and calculates cell model average
    :param base: Path to base image
    :param fluor: Path to fluor image
    :param basetype: One of Membrane, BF or Phase
    :return:
    """
    if not os.path.exists(fluor):
        return 0

    par = ParametersManager()

    if basetype == "Phase":
        par.imageloaderparams.mask_fill_holes = True
    elif basetype == "BF":
        par.imageloaderparams.mask_algorithm = 'StarDist_BF'
    elif basetype == "Membrane":
        par.imageloaderparams.mask_algorithm = 'StarDist'
    else:
        raise ValueError(f"Only Phase, BF or Membrane are valid base image types. You provided {basetype}")

    # Load images and compute mask
    imm = ImageManager()
    imm.load_base_image(base, par.imageloaderparams)
    imm.compute_base_mask(par.imageloaderparams)
    imm.compute_mask(par.imageloaderparams)

    if os.path.exists(membrane):
        imm.load_fluor_image(membrane, par.imageloaderparams)
        imm.load_option_image(fluor, par.imageloaderparams)
    else:
        imm.load_fluor_image(fluor, par.imageloaderparams)

    seg = SegmentsManager()
    seg.compute_segments(par.imageprocessingparams, imm)

    seg.save_labels(filename=os.path.join(os.path.dirname(base), 'Labels'))

    cel = CellManager(par)
    cel.compute_cells(par, imm, seg)
    cel.process_cells(par.cellprocessingparams, imm)
    cel.overlay_cells(imm)

    if os.path.exists(membrane):
        print("BLALBLAL")
        # Classify cells
        ccc = CellCycleClassifier()
        ccc.classify_cells(imm, cel, "Epifluorescence", False)

        # Cell alignment
        cellalligner = CellAligner(cel, fluorescence="Optional")
        cellalligner.align_cells()

    else:
        # Cell alignment
        cellalligner = CellAligner(cel)
        cellalligner.align_cells()

    # Cell model
    savepath = fluor[:-4] + "_CA_"
    cellmodeler = CellModeler(cel)

    if os.path.exists(fluor.replace('.tif', '.xml')) and os.path.exists(membrane):
        cellmodeler.create_cell_model_from_TM(fluor.replace('.tif', '.xml'), savepath, imm.align_values, True)
    elif os.path.exists(fluor.replace('.tif', '.xml')) and not os.path.exists(membrane):
        cellmodeler.create_cell_model_from_TM(fluor.replace('.tif', '.xml'), savepath, imm.align_values)
    else:
        cellmodeler.create_cell_model_from_ehooke(savepath)

    return 1


def run_batch_workflow(folder, fluor_names, base_name, base_type, membrane):
    for replicate in os.listdir(folder):
        replicate_fullname = os.path.join(folder, replicate)
        if os.path.isdir(replicate_fullname):
            base_path = os.path.join(replicate_fullname, base_name + '.tif')
            for fname in fluor_names:
                fluor_path = os.path.join(replicate_fullname, fname + '.tif')
                membrane_path = os.path.join(replicate_fullname, membrane + '.tif')
                run_CellAverager(base_path, fluor_path, base_type, membrane_path)


def join_replicates(folder, fluor_names, phase):

    all_fluor1 = []
    N_fluor1 = []
    all_fluor1_TM = []
    N_fluor1_TM = []

    all_fluor2 = []
    N_fluor2 = []
    all_fluor2_TM = []
    N_fluor2_TM = []

    for replicate in os.listdir(folder):
        replicate_fullname = os.path.join(folder, replicate)

        if os.path.isdir(replicate_fullname):
            fluor1name = os.path.join(replicate_fullname, fluor_names[0]+f"_CA_model_greater_1_phase{phase}"+".tif")
            fluor1name_tm = os.path.join(replicate_fullname, fluor_names[0]+f"_CA_model_greater_1_phase{phase}_TM_SPOTS"+".tif")
            fluor2name = os.path.join(replicate_fullname, fluor_names[1]+f"_CA_model_greater_1_phase{phase}"+".tif")
            fluor2name_tm = os.path.join(replicate_fullname, fluor_names[1]+f"_CA_model_greater_1_phase{phase}_TM_SPOTS"+".tif")

            if os.path.exists(fluor1name):
                with tifffile.TiffFile(fluor1name) as tif:
                    all_fluor1.append(tif.asarray())
                    N_fluor1.append(int(tif.pages[0].tags["image_description"].value))
                with tifffile.TiffFile(fluor1name_tm) as tif:
                    all_fluor1_TM.append(tif.asarray())
                    N_fluor1_TM.append(int(tif.pages[0].tags["image_description"].value))

            if os.path.exists(fluor2name):
                with tifffile.TiffFile(fluor2name) as tif:
                    all_fluor2.append(tif.asarray())
                    N_fluor2.append(int(tif.pages[0].tags["image_description"].value))
                with tifffile.TiffFile(fluor2name_tm) as tif:
                    all_fluor2_TM.append(tif.asarray())
                    N_fluor2_TM.append(int(tif.pages[0].tags["image_description"].value))

    cumsum_fluor1 = np.zeros((30, 30))
    for idx, N in enumerate(N_fluor1):
        cumsum_fluor1 += resize(all_fluor1[idx] * N, (30, 30))
    cumsum_fluor2 = np.zeros((30, 30))
    for idx, N in enumerate(N_fluor2):
        cumsum_fluor2 += resize(all_fluor2[idx] * N, (30, 30))

    cumsum_fluor1_TM = np.zeros((30, 30))
    for idx, N in enumerate(N_fluor1_TM):
        cumsum_fluor1_TM += resize(all_fluor1_TM[idx] * N, (30, 30))
    cumsum_fluor2_TM = np.zeros((30, 30))
    for idx, N in enumerate(N_fluor2_TM):
        cumsum_fluor2_TM += resize(all_fluor2_TM[idx] * N, (30, 30))

    cellmodeler = CellModeler(None)
    colormap = matplotlib.cm.get_cmap("coolwarm")

    if len(N_fluor1) > 0:
        imsave(folder + os.sep + f"model_{fluor_names[0]}_{phase}.tif", cumsum_fluor1 / sum(N_fluor1), plugin="tifffile")
        color_img1 = np.zeros((30, 30, 3))
        color_fluor1 = cellmodeler.assign_color(cumsum_fluor1 / sum(N_fluor1), color_img1, colormap)
        imsave(folder + os.sep + f"color_{fluor_names[0]}_n{sum(N_fluor1)}_{phase}.png", color_fluor1)

    if len(N_fluor1_TM) > 0:
        imsave(folder + os.sep + f"model_{fluor_names[0]}_TM__{phase}.tif", cumsum_fluor1_TM / sum(N_fluor1_TM), plugin="tifffile")
        color_img1_TM = np.zeros((30, 30, 3))
        color_fluor1_TM = cellmodeler.assign_color(cumsum_fluor1_TM / sum(N_fluor1_TM), color_img1_TM, colormap)
        imsave(folder + os.sep + f"color_{fluor_names[0]}_TM_n{sum(N_fluor1_TM)}_{phase}.png", color_fluor1_TM)

    if len(N_fluor2) > 0:
        imsave(folder + os.sep + f"model_{fluor_names[1]}_{phase}.tif", cumsum_fluor2 / sum(N_fluor2), plugin="tifffile")
        color_img2 = np.zeros((30, 30, 3))
        color_fluor2 = cellmodeler.assign_color(cumsum_fluor2 / sum(N_fluor2), color_img2, colormap)
        imsave(folder + os.sep + f"color_{fluor_names[1]}_n{sum(N_fluor2)}_{phase}.png", color_fluor2)

    if len(N_fluor2_TM) > 0:
        imsave(folder + os.sep + f"model_{fluor_names[1]}_TM__{phase}.tif", cumsum_fluor2_TM / sum(N_fluor2_TM), plugin="tifffile")
        color_img2_TM = np.zeros((30, 30, 3))
        color_fluor2_TM = cellmodeler.assign_color(cumsum_fluor2_TM / sum(N_fluor2_TM), color_img2_TM, colormap)
        imsave(folder + os.sep + f"color_{fluor_names[1]}_TM_n{sum(N_fluor2_TM)}_{phase}.png", color_fluor2_TM)


if __name__ == '__main__':

    gui = MainGUI()
    gui.mainloop()

    parameters = gui.pars

    ehooke = parameters['ehookefolder']
    rootfolder = parameters['rootfolder']
    fluornames = [parameters['fluor1'], parameters['fluor2']]
    basename = parameters['base']
    basetype = parameters['basetype']
    membrane = parameters['membrane']

    sys.path.append(ehooke)

    from images import ImageManager
    from parameters import MaskParameters, RegionParameters, CellParameters, ParametersManager
    from segments import SegmentsManager
    from cells import CellManager
    from cellcycleclassifier import CellCycleClassifier

    for condition in os.listdir(rootfolder):

        path = os.path.join(rootfolder, condition)

        if not os.path.isdir(path):
            continue

        print(condition)
        run_batch_workflow(path, fluornames, basename, basetype, membrane)
        join_replicates(path, fluornames, (0,))
