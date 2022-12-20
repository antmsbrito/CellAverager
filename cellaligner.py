import os
import math
import numpy as np
from skimage.io import imsave
from skimage.transform import rotate
from skimage.morphology import binary_erosion
from sklearn.decomposition import PCA
from tkinter import filedialog as fd


class CellAligner:
    """
    This class is in charge of taking a cellmanager object from ehooke and align all cells
    returns a cell manager object with cell objects with 2 new attributes
        - aligned cell mask
        - aligned fluor mask
    """

    def __init__(self, cellmanagerobj):
        self.cellmanager = cellmanagerobj

    def calculate_rotation_angle(self, ids):
        binary = self.cellmanager.cells[ids].fluor * self.cellmanager.cells[ids].cell_mask
        outline = self.calculate_cell_outline(binary)
        major_axis = self.calculate_major_axis(outline)
        return self.calculate_axis_angle(major_axis)

    def align_cells(self):

        for key in self.cellmanager.cells:
            angle = self.calculate_rotation_angle(key)
            self.cellmanager.cells[key].angle_of_rotation = angle
            self.cellmanager.cells[key].aligned_cell_mask = rotate(self.cellmanager.cells[key].cell_mask, angle)
            self.cellmanager.cells[key].aligned_fluor_mask = rotate(
                self.cellmanager.cells[key].fluor * self.cellmanager.cells[key].cell_mask, angle)

    def save_aligned_cells(self, path=None):
        if path is None:
            path = fd.askdirectory(title="SAVE ALIGNED CELLS")

        for key in self.cellmanager.cells:
            imsave(path + os.sep + key + ".tif", self.cellmanager.cells[key].aligned_fluor_mask)

    @staticmethod
    def calculate_cell_outline(binary):
        outline = binary * (1 - binary_erosion(binary))

        return outline

    @staticmethod
    def calculate_major_axis(outline):
        x, y = np.nonzero(outline)
        x = [[val] for val in x]
        y = [[val] for val in y]
        coords = np.concatenate((x, y), axis=1)

        pca = PCA(n_components=1)
        pca.fit(coords)

        pos_x, pos_y = pca.mean_
        eigenvector_x, eigenvector_y = pca.components_[0]
        eigenval = pca.explained_variance_[0]

        return [[pos_x - eigenvector_x * eigenval, pos_y - eigenvector_y * eigenval],
                [pos_x + eigenvector_x * eigenval, pos_y + eigenvector_y * eigenval]]

    @staticmethod
    def calculate_axis_angle(major_axis):
        # TODO refactor, atan2 picks correct quadrant
        x0, y0 = major_axis[0]
        x1, y1 = major_axis[1]

        if x0 - x1 == 0:
            angle = 0.0

        elif y0 - y1 == 0:
            angle = 90.0

        else:
            if y1 > y0:
                if x1 > x0:
                    direction = -1
                    opposite = x1 - x0
                    adjacent = y1 - y0
                else:
                    direction = 1
                    opposite = x0 - x1
                    adjacent = y1 - y0

            elif y0 > y1:
                if x1 > x0:
                    direction = 1
                    opposite = x1 - x0
                    adjacent = y0 - y1
                else:
                    direction = -1
                    opposite = x0 - x1
                    adjacent = y0 - y1

            angle = math.degrees(math.atan(opposite / adjacent)) * direction

        if angle != 0:
            angle = 90.0 - angle
        else:
            angle = 90

        return angle