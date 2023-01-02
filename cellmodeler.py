import os
import numpy as np
import matplotlib as mpl
from tkinter import filedialog as fd
import xml.etree.ElementTree as ET
from skimage.io import imsave
from skimage.transform import resize
from skimage.filters import threshold_isodata
from skimage.color import gray2rgb

from skimage.draw import circle
from skimage.transform import rotate

from cellspots import Spots

from matplotlib import pyplot as plt


class CellModeler:
    """
    This class is in charge of taking a cellmanager object from ehooke THAT HAS BEEN pre aligned
    returns

    """

    def __init__(self, cellmanager):

        self.cellmanager = cellmanager

        self.cell_model = None
        self.color_model = None
        self.mean_x = 0
        self.mean_y = 0
        self.number_of_cells = None

    def resize_cells(self, cells):
        tmp = []
        for cell in cells:
            tmp.append(resize(cell.aligned_fluor_mask, (self.mean_x, self.mean_y)))

        return np.array(tmp)

    def create_cell_average(self, cells):
        model_cell = np.zeros((self.mean_x, self.mean_y))
        for cell in np.array(cells):
            model_cell += cell

        model_cell /= float(len(cells))

        return model_cell

    def create_cell_model_from_ehooke(self, savepath):

        self.mean_x = int(
            np.median([self.cellmanager.cells[key].aligned_fluor_mask.shape[0] for key in self.cellmanager.cells]))
        self.mean_y = int(
            np.median([self.cellmanager.cells[key].aligned_fluor_mask.shape[1] for key in self.cellmanager.cells]))

        self.number_of_cells = len(self.cellmanager.cells)
        selected_cells = [self.cellmanager.cells[key] for key in self.cellmanager.cells]
        selected_cells = self.resize_cells(selected_cells)
        self.cell_model = self.create_cell_average(selected_cells)

        self.save_cell_model(savepath + f"model")
        self.model2color(savepath + f"color")

    def create_cell_model_from_TM(self, xmlfile, savepath, align):

        spots = self.read_xml_spots(xmlfile)

        for key in self.cellmanager.cells:
            box = self.cellmanager.cells[key].box
            self.cellmanager.cells[key].spots, self.cellmanager.cells[key].spots_coords = spots.filterbox(box, align)

        self.create_cell_model(1, 'greater', savepath)

    def create_cell_model(self, spotnumber, operation, savepath):

        if operation == 'equal':
            selected_cells = [self.cellmanager.cells[key] for key in self.cellmanager.cells if
                              self.cellmanager.cells[key].spots == spotnumber]
        elif operation == 'greater':
            selected_cells = [self.cellmanager.cells[key] for key in self.cellmanager.cells if
                              self.cellmanager.cells[key].spots >= spotnumber]
        elif operation == 'smaller':
            selected_cells = [self.cellmanager.cells[key] for key in self.cellmanager.cells if
                              self.cellmanager.cells[key].spots <= spotnumber]
        else:
            raise (ValueError(f"Unrecognized operation: {operation}"))

        if len(selected_cells) == 0:
            print(f"No cells with {spotnumber} spots ", operation)
            return 0

        self.mean_x = int(np.median([s.aligned_fluor_mask.shape[0] for s in selected_cells]))
        self.mean_y = int(np.median([s.aligned_fluor_mask.shape[1] for s in selected_cells]))
        self.number_of_cells = len(selected_cells)

        new_cells = []
        for cell in selected_cells:
            blank = np.zeros(cell.aligned_fluor_mask.shape)
            for c in cell.spots_coords:
                xc, yc = c
                rr, cc = circle(yc, xc, 1, shape=blank.shape)
                blank[rr, cc] = 1
            blank = rotate(blank, cell.angle_of_rotation)
            new_cells.append(blank)

        selected_cells = self.resize_cells(selected_cells)
        self.cell_model = self.create_cell_average(selected_cells)
        self.save_cell_model(savepath + f"model_{operation}_{spotnumber}")
        self.model2color(savepath + f"color_{operation}_{spotnumber}", selected_cells)

        new_cells = np.array([resize(c, (self.mean_x, self.mean_y)) for c in new_cells])
        self.cell_model = self.create_cell_average(new_cells)
        self.save_cell_model(savepath + f"model_{operation}_{spotnumber}_TM_SPOTS")
        self.model2color(savepath + f"color_{operation}_{spotnumber}_TM_SPOTS", new_cells)

    def save_cell_model(self, path=None):
        if path is None:
            path = fd.asksaveasfilename()

        imsave(path + ".tif", self.cell_model, plugin="tifffile", imagej=False, description=str(self.number_of_cells))

    def model2color(self, savepath, celllist):

        mask = self.cell_model > threshold_isodata(self.cell_model)
        filtered = self.cell_model * mask

        colormap = mpl.cm.get_cmap("coolwarm")
        color_img = np.zeros(np.shape(gray2rgb(filtered)))

        self.color_model = self.assign_color(filtered, color_img, colormap)

        imsave(savepath + f"_n{self.number_of_cells}.png", self.color_model)

    @staticmethod
    def read_xml_spots(xmlfile):
        # Read xml file into memory
        root = ET.parse(xmlfile).getroot()
        model_child = root[0]  # This is trackmate model object
        allspots = model_child[1]  # This is the spot xml node
        quality = float(root[1][4][0].attrib['value'])

        return Spots(allspots, quality)

    @staticmethod
    def assign_color(modelmasked, outimage, cmap):
        norm = mpl.colors.Normalize(vmin=np.amin(modelmasked[np.nonzero(modelmasked)]), vmax=np.amax(modelmasked))
        for i in range(modelmasked.shape[0]):
            for ii in range(modelmasked.shape[1]):
                px = modelmasked[i, ii]
                if np.abs(px) < 1e-3:
                    outimage[i, ii] = (0, 0, 0)
                else:
                    rgba = cmap(norm(px))
                    outimage[i, ii] = (rgba[0], rgba[1], rgba[2])

        return outimage
