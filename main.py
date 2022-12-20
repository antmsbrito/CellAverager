import os
import sys

sys.path.append(r'C:\Users\António\Documents\GitHub\eHooke_1.0')

import numpy as np
import matplotlib
from skimage.transform import resize
from skimage.external import tifffile
from skimage.io import imsave
from skimage.color import gray2rgb

from images import ImageManager
from parameters import MaskParameters, RegionParameters, CellParameters, ParametersManager
from segments import SegmentsManager
from cells import CellManager

from cellmodeler import CellModeler
from cellaligner import CellAligner


def runCA(phase, fluor, savepath):
    """
    Given a set of images run ehooke analysis
    """
    if not os.path.exists(fluor):
        return 0

    # Set params
    par = ParametersManager()

    par.imageloaderparams.mask_fill_holes = True

    # Load images
    imm = ImageManager()
    imm.load_base_image(phase, par.imageloaderparams)
    imm.compute_base_mask(par.imageloaderparams)
    imm.compute_mask(par.imageloaderparams)
    imm.load_fluor_image(fluor, par.imageloaderparams)

    seg = SegmentsManager()
    seg.compute_segments(par.imageprocessingparams, imm)

    seg.save_labels(filename=os.path.join(os.path.dirname(phase), 'Labels'))

    # TODO CellManager and process cells needs parameters manager instance
    # All other functions do not and instead take parameters individually
    # Fix inconsistency on ehooke source code, either take parmanager or individual pars
    cel = CellManager(par)
    cel.compute_cells(par, imm, seg)
    cel.process_cells(par.cellprocessingparams, imm)
    cel.overlay_cells(imm)

    cellalligner = CellAligner(cel)
    cellalligner.align_cells()

    cellmodeler = CellModeler(cel)
    cellmodeler.create_cell_model_from_TM(fluor.replace('.tif', '.xml'), savepath, imm.align_values)

    return len(cel.cells)


def run_workflow(folder):

    ## Terminus
    nt = runCA(folder + r"\Phase.tif",
               folder + r"\CFP.tif",
               folder + r"\CFP_hm")

    ## Origin
    no = runCA(folder + r"\Phase.tif",
               folder + r"\YFP.tif",
               folder + r"\YFP_hm")

    return nt


def run_batch_workflow(folder):
    total_cells = 0
    for replicate in os.listdir(folder):
        replicate_fullname = os.path.join(folder, replicate)
        if os.path.isdir(replicate_fullname):
            n = run_workflow(replicate_fullname)
            total_cells += n

    print(folder, total_cells)


def join_replicates(folder, spot_condition):
    all_cfp = []
    all_yfp = []
    N_cfp = []
    N_yfp = []

    for replicate in os.listdir(folder):
        replicate_fullname = os.path.join(folder, replicate)

        if os.path.isdir(replicate_fullname):
            if os.path.exists(os.path.join(replicate_fullname, "CFP_hm" + spot_condition + ".tif")):
                with tifffile.TiffFile(replicate_fullname + os.sep + "CFP_hm" + spot_condition + ".tif") as tif:
                    all_cfp.append(tif.asarray())
                    N_cfp.append(int(tif.pages[0].tags["image_description"].value))

            if os.path.exists(os.path.join(replicate_fullname, "YFP_hm" + spot_condition + ".tif")):
                with tifffile.TiffFile(replicate_fullname + os.sep + "YFP_hm" + spot_condition + ".tif") as tif:
                    all_yfp.append(tif.asarray())
                    N_yfp.append(int(tif.pages[0].tags["image_description"].value))

    cumsum_cfp = np.zeros((30, 30))
    for idx, N in enumerate(N_cfp):
        cumsum_cfp += resize(all_cfp[idx] * N, (30, 30))

    cumsum_yfp = np.zeros((30, 30))
    for idx, N in enumerate(N_yfp):
        cumsum_yfp += resize(all_yfp[idx] * N, (30, 30))

    cellmodeler = CellModeler(None)
    colormap = matplotlib.cm.get_cmap("coolwarm")

    if len(N_cfp) > 0:
        imsave(folder + os.sep + f"model_cfp_{spot_condition}.tif", cumsum_cfp / sum(N_cfp), plugin="tifffile")
        color_img = np.zeros((30, 30, 3))
        color_cfp = cellmodeler.assign_color(cumsum_cfp / sum(N_cfp), color_img, colormap)
        imsave(folder + os.sep + f"color_cfp_{spot_condition}_n{sum(N_cfp)}.png", color_cfp)

    if len(N_yfp) > 0:
        imsave(folder + os.sep + f"model_yfp_{spot_condition}.tif", cumsum_yfp / sum(N_yfp), plugin="tifffile")
        color_img = np.zeros((30, 30, 3))
        color_yfp = cellmodeler.assign_color(cumsum_yfp / sum(N_yfp), color_img, colormap)

        imsave(folder + os.sep + f"color_yfp_{spot_condition}_n{sum(N_yfp)}.png", color_yfp)


if __name__ == '__main__':

    for condition in os.listdir(r"C:\Users\António\Documents\30-11-2022_heatmaps"):

        path = os.path.join(r"C:\Users\António\Documents\30-11-2022_heatmaps", condition)
        print(condition)
        run_batch_workflow(path)
        join_replicates(path, "model_greater_1_TM_SPOTS")
        join_replicates(path, "model_greater_1")
