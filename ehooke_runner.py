import sys

sys.path.append(r'C:\Users\António\Desktop\eHooke_1.0-master')

from images import ImageManager
from parameters import MaskParameters, RegionParameters, CellParameters, ParametersManager
from segments import SegmentsManager
from cells import CellManager

from cellmodeler import CellModeler
from cellaligner import CellProcessor
from process_cell_model import model2color

import os
from PIL import Image
from skimage.io import imsave
from skimage.morphology import binary_dilation
from skimage.exposure import rescale_intensity
from skimage.external import tifffile
from skimage.transform import resize
import numpy as np


def run_ehooke(phase, fluor):
    """
    Given a set of images run ehooke analysis
    """

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

    seg.save_labels(filename=os.path.join(os.path.dirname(phase), 'labels'))

    #TODO CellManager and process cells needs parameters manager instance
    # All other functions do not and instead take parameters individually
    # Fix inconsistency on ehooke source code, either take parmanager or individual pars
    cel = CellManager(par)
    cel.compute_cells(par, imm, seg)
    cel.process_cells(par.cellprocessingparams, imm)
    cel.overlay_cells(imm)

    return imm, seg, cel


def run_mean_cell(im, ce, savepath):
    crops_original = {}
    mask_dict = {}
    for key in ce.cells:
        cell = ce.cells[key]
        mask = cell.cell_mask

        count = 0
        while count < -1:
            count += 1
            mask = binary_dilation(mask)

        # ORIGINAL FLUOR
        x0, y0, x1, y1 = cell.box
        ori_fluor = im.original_fluor_image[x0:x1 + 1, y0:y1 + 1]
        crops_original[key] = ori_fluor * mask

        mask_dict[key] = mask

    cellalligner = CellProcessor(crops_original, mask_dict)
    cellalligner.align_cells()
    # cellalligner.save_aligned_cells(r"original")
    cellmodeler = CellModeler()
    cellmodeler.create_cell_model_from_ehooke(cellalligner.aligned_cells)
    cellmodeler.save_cell_model(savepath)


def run_workflow(folder):
    ## ORIGIN
    im, se, ce = run_ehooke(folder + r'\Phase.tif',
                            folder + r'\YFP.tif')
    run_mean_cell(im, ce, folder + r'\model_origin')
    model2color(folder + r'\model_origin.tif',
                folder + r'\color_origin.png')
    ## TERMINUS
    im, se, ce = run_ehooke(folder + r'\Phase.tif',
                            folder + r'\CFP.tif')
    run_mean_cell(im, ce, folder + r'\model_terminus')
    model2color(folder + r'\model_terminus.tif',
                folder + r'\color_terminus.png')


def run_batch_workflow(folder):
    for replicate in os.listdir(folder):
        replicate_fullname = os.path.join(folder, replicate)
        if os.path.isdir(replicate_fullname):
            run_workflow(replicate_fullname)


def join_replicates(folder):
    all_origins = []
    all_terminus = []
    N_origins = []
    N_terminus = []

    for replicate in os.listdir(folder):
        replicate_fullname = os.path.join(folder, replicate)
        if os.path.isdir(replicate_fullname):
            with tifffile.TiffFile(replicate_fullname + os.sep + "model_origin.tif") as tif:
                all_origins.append(tif.asarray())
                N_origins.append(int(tif.pages[0].tags["image_description"].value))

            with tifffile.TiffFile(replicate_fullname + os.sep + "model_terminus.tif") as tif:
                all_terminus.append(tif.asarray())
                N_terminus.append(int(tif.pages[0].tags["image_description"].value))

            assert N_origins[-1] == N_terminus[-1]

    # origins
    cumsum_ori = np.zeros((30, 30))
    cumsum_ter = np.zeros((30, 30))
    for idx, N in enumerate(N_origins):
        cumsum_ter += resize(all_terminus[idx] * N, (30, 30))
        cumsum_ori += resize(all_origins[idx] * N, (30, 30))

    imsave(folder + os.sep + "model_terminus.tif", cumsum_ter / sum(N_origins), plugin="tifffile")
    model2color(folder + r'\model_terminus.tif', folder + r'\color_terminus.png')
    imsave(folder + os.sep + "model_origin.tif", cumsum_ori / sum(N_origins), plugin="tifffile")
    model2color(folder + r'\model_origin.tif', folder + r'\color_origin.png')

    print(folder, sum(N_origins))


if __name__ == '__main__':
    """
    run_workflow(r"imgs_test\froswt")
    run_workflow(r"imgs_test\fros_d0")
    run_workflow(r"imgs_test\fros_dscp")
    run_workflow(r"imgs_test\fros_dspo0j")
  

    run_batch_workflow(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dSAUSA300_0383")
    run_batch_workflow(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dsbcE")
    run_batch_workflow(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dscpAB")
    run_batch_workflow(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dspo0J")
    run_batch_workflow(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\JE2")
    """

    join_replicates(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dSAUSA300_0383")
    join_replicates(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dsbcE")
    join_replicates(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dscpAB")
    join_replicates(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\dspo0J")
    join_replicates(r"E:\SC_spots\FROS mutants\21-7-2022_JE2-FROS-mutants\JE2")
