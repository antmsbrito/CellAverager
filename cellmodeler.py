import xml.etree.ElementTree as ET

import numpy as np
import matplotlib as mpl
from skimage.draw import circle
from skimage.transform import resize, rotate

from cellspots import Spots

class CellModeler:
    """
    Class that manages and creates cell models based on eHooke classes THAT HAVE BEEN PRE ALIGNED
    There is support for an optional trackmate xml with spot detection to filter for spots
    By default always tries to compute models with cells with more than 1 spot 
    By default always checks for cell cycle
    """
    def __init__(self, cellmanager, imagemanager, xmlfile:str|None=None) -> None:
        
        # Store ehooke classes 
        self.cellmanager = cellmanager
        self.imagemanager = imagemanager

        # Check if xmlfile exists and if it does store the number of spots for each cell
        if xmlfile:
            self.spots = self.read_xml_spots(xmlfile)
            for key in self.cellmanager.cells:
                box = self.cellmanager.cells[key].box
                self.cellmanager.cells[key].spots, self.cellmanager.cells[key].spots_coords = self.spots.filterbox(box, self.imagemanager.align_values)
        else:
            self.spots = None

            
    def select_cells(self, minspots:int|np.inf=1, maxspots:int|np.inf=np.inf, cellcycle:tuple=(0,1,2,3)) -> list:
        """
        This function builds a selection of cells based on spot number and desired cell cycle phase
        In the case of no spot detection it ignores the min and max spot given
        In the case of no cell cycle phase, all cells are phase 0
        """
        
        if self.spots:
            selection = [self.cellmanager.cells[key] for key in self.cellmanager.cells if
                        (maxspots>self.cellmanager.cells[key].spots>minspots) and 
                        (self.cellmanager.cells[key].stats["Cell Cycle Phase"] in cellcycle)]
        else:
            selection = [self.cellmanager.cells[key] for key in self.cellmanager.cells if 
                        (self.cellmanager.cells[key].stats["Cell Cycle Phase"] in cellcycle)]

        return selection

    def build_spot_model(self, selected_cells:list) -> np.ndarray:
        """
        This functions builds a cell model based on spot localization within each cell
        Returns the empty array if there is not spot detection xml
        """

        if not self.spots:
            return np.array([])

        spot_cells = []
        for cell in selected_cells:
            blank = np.zeros(cell.aligned_fluor_mask.shape)
            for c in cell.spots_coords:
                xc, yc = c
                rr, cc = circle(yc, xc, 1, shape=blank.shape)
                blank[rr, cc] = 1
            blank = rotate(blank, cell.angle_of_rotation)
            spot_cells.append(blank)

        x_size = int(np.median([s.shape[0] for s in spot_cells]))
        y_size = int(np.median([s.shape[1] for s in spot_cells]))

        resized_cells = self.resize_arr(spot_cells, x_size, y_size)
        cell_model = self.create_average(resized_cells)

        return cell_model

    def build_average_model(self, selected_cells:list) -> np.ndarray:
        """
        This function builds a cell model based on the average of the aligned fluorescence 
        """

        aligned_fluor = [selected_cells[key].aligned_fluor_mask for key in selected_cells]

        x_size = int(np.median([s.shape[0] for s in aligned_fluor]))
        y_size = int(np.median([s.shape[1] for s in aligned_fluor]))

        resized_cells = self.resize_arr(aligned_fluor, x_size, y_size)
        cell_model = self.create_average(resized_cells)

        return cell_model

    @staticmethod
    def resize_arr(imgarr:list, xsize:int, ysize:int)->np.ndarray:
        resized = np.zeros(xsize, ysize, len(imgarr))
        for idx, img in enumerate(imgarr):
            resized[:,:,idx] = resize(img, (xsize, ysize))
        return resized

    @staticmethod
    def create_average(imgarr:np.ndarray)->np.ndarray:
        return np.average(imgarr, axis=2)
        
    @staticmethod
    def read_xml_spots(xmlfile:str)->Spots:
        # Read xml file into memory
        root = ET.parse(xmlfile).getroot()
        model_child = root[0]  # This is trackmate model object
        allspots = model_child[1]  # This is the spot xml node
        quality = float(root[1][4][0].attrib['value'])

        return Spots(allspots, quality)