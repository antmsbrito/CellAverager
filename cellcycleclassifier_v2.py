import numpy as np
from skimage.util import img_as_float
from skimage.exposure import rescale_intensity
from skimage.transform import resize as skresize
from keras.models import load_model

from scipy.signal import fftconvolve
from scipy.ndimage import center_of_mass
from skimage.transform import EuclideanTransform, warp
from skimage.io import imread

# force classification to happen on CPU to avoid CUDA problems
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

class CellCycleClassifier(object):

    def __init__(self, membrane, dnapath, imagemanager, cellmanager):

        self.model = load_model("cellcycle_cnn_model")
        
        self.membrane = membrane
        self.dna = self.load_option_image(dnapath, imagemanager.mask)
        self.cell_manager = cellmanager

    def preprocess_image(self, image, microscope):

        if microscope == "Epifluorescence":
            max_dim = 50

        elif microscope == "SIM":
            max_dim = 100

        h, w = image.shape

        max_h, max_w = max_dim, max_dim

        lines_to_add = max_h - h
        columns_to_add = max_w - w

        if lines_to_add > 0:
            if lines_to_add%2 == 0:
                new_line = np.zeros((int(lines_to_add/2), w))
                image = np.concatenate((new_line, image, new_line), axis=0)
            else:
                new_line_top = np.zeros((int(lines_to_add/2)+1, w))
                new_line_bot = np.zeros((int(lines_to_add/2), w))
                image = np.concatenate((new_line_top, image, new_line_bot), axis=0)

        elif lines_to_add < 0:
            if (lines_to_add * -1)%2 == 0:
                cutsize = int((lines_to_add * -1)/2)
                image = image[cutsize:h-cutsize, :]
            else:
                cutsize = int((lines_to_add * -1)/2)
                image = image[cutsize:h-cutsize-1, :]


        if columns_to_add > 0:
            if columns_to_add%2 == 0:
                columns_to_add = np.zeros((max_dim, int(columns_to_add/2)))
                image = np.concatenate((columns_to_add, image, columns_to_add), axis=1)
            else:
                columns_to_add_left = np.zeros((max_dim, int(columns_to_add/2)+1))
                columns_to_add_right = np.zeros((max_dim, int(columns_to_add/2)))
                image = np.concatenate((columns_to_add_left, image, columns_to_add_right), axis=1)

        elif columns_to_add < 0:
            if (columns_to_add * -1)%2 == 0:
                cutsize = int((columns_to_add * -1)/2)
                image = image[:, cutsize:w-cutsize]
            else:
                cutsize = int((columns_to_add * -1)/2)
                image = image[:, cutsize:w-cutsize-1]

        image = img_as_float(image)
        image = image.reshape(max_dim, max_dim, 1)

        return image

    def classify_cell(self, fluor, optional, microscope):

        fluor_img = skresize(self.preprocess_image(fluor, microscope),
                             (100, 100),
                             order=0,
                             preserve_range=True,
                             anti_aliasing=False,
                             anti_aliasing_sigma=None)

        optional_img = skresize(self.preprocess_image(optional, microscope),
                                (100, 100),
                                order=0,
                                preserve_range=True,
                                anti_aliasing=False,
                                anti_aliasing_sigma=None)

        pred = self.model.predict_classes(np.concatenate((fluor_img, optional_img), axis=1).reshape(-1, 100, 200, 1))

        return pred[0] + 1

    def classify_cells(self, microscope):
        
        for k in self.cell_manager.cells.keys():
            cell = self.cell_manager.cells[k]

            x0, y0, x1, y1 = cell.box

            cell_fluor = rescale_intensity(self.membrane[x0:x1 + 1, y0:y1 + 1] * cell.cell_mask)
            cell_optional = rescale_intensity(self.dna[x0:x1 + 1, y0:y1 + 1] * cell.cell_mask)

            self.cell_manager.cells[k].stats["Cell Cycle Phase"] = self.classify_cell(cell_fluor, cell_optional, microscope)
            
            
    def load_option_image(self, filename, mask):
        """Loads an option image that can be used to look for the septum
        and to help classify the cell cycle phases. No fluorescence is measured
        on this image"""

        inverted_mask = 1 - mask

        optional_image = imread(filename)

        if len(optional_image.shape) > 2:
            optional_image = color.rgb2gray(optional_image)

        optional_image = img_as_float(optional_image)

        best = (0, 0)
        
        # Alignment is done by taking the maximum of the correlation
        # between phase and fluorescence
        corr = fftconvolve(inverted_mask, optional_image[::-1, ::-1])
        deviation = np.unravel_index(np.argmax(corr), corr.shape)
        cm = center_of_mass(np.ones(corr.shape))
        best = np.subtract(deviation, cm)

        dx, dy = best
        matrix = EuclideanTransform(rotation=0, translation=(dx, dy))
        return warp(optional_image, matrix.inverse, preserve_range=True)