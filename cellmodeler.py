import xml.etree.ElementTree as ET
import math
import numpy as np
import matplotlib as mpl
from skimage.draw import circle
from skimage.transform import resize, rotate

#from cellspots import Spots
from cellspotsNewXML import Spots as SpotsV2

class CellModeler:
    """
    Class that manages and creates cell models based on eHooke classes THAT HAVE BEEN PRE ALIGNED
    There is support for an optional trackmate xml with spot detection to filter for spots
    By default always tries to compute models with cells with more than 1 spot 
    By default always checks for cell cycle
    """
    def __init__(self, cellmanager, imagemanager, xmlfile=None) -> None:
        
        # Store ehooke classes 
        self.cellmanager = cellmanager
        self.imagemanager = imagemanager
        
        self.total_cells = len(self.cellmanager.cells)

        # Check if xmlfile exists and if it does store the number of spots for each cell
        if xmlfile:
            self.spots = self.read_xml_spots(xmlfile)
            for key in self.cellmanager.cells:
                if self.cellmanager.cells[key].axis_ratio:
                    box = self.cellmanager.cells[key].box
                    #self.cellmanager.cells[key].spots, self.cellmanager.cells[key].spots_coords = self.spots.filterbox(box, self.imagemanager.align_values)
                    # TODO stardistlabelsmight not exist
                    self.cellmanager.cells[key].spots, self.cellmanager.cells[key].spots_coords = self.spots.filterlabel(box,self.cellmanager.cells[key].label,self.imagemanager.stardist_labels)
                else:
                    self.cellmanager.cells[key].spots = -1
        else:
            self.spots = None

            
    def select_cells(self, minspots=1, maxspots=np.inf, cellcycle:tuple=(0,1,2,3)) -> list:
        """
        This function builds a selection of cells based on spot number and desired cell cycle phase
        In the case of no spot detection it ignores the min and max spot given
        In the case of no cell cycle phase, all cells are phase 0
        """
        
        if self.spots:
            selection = [self.cellmanager.cells[key] for key in self.cellmanager.cells if
                        (maxspots>self.cellmanager.cells[key].spots>=minspots) and 
                        (self.cellmanager.cells[key].stats["Cell Cycle Phase"] in cellcycle)]
            n_spots = sum([s.spots for s in selection])

        else:
            selection = [self.cellmanager.cells[key] for key in self.cellmanager.cells if 
                        (self.cellmanager.cells[key].stats["Cell Cycle Phase"] in cellcycle)]
            n_spots = 0
        
        return selection, n_spots

    def build_spot_model(self, selected_cells:list) -> np.ndarray:
        """
        This functions builds a cell model based on spot localization within each cell
        Returns the empty array if there is not spot detection xml
        """

        if not self.spots:
            return np.array([])
        if len(selected_cells) == 0:
            return np.array([])

        spot_cells = []
        for cell in selected_cells:
            blank = np.zeros([cell.aligned_fluor_mask.shape[0]-8,cell.aligned_fluor_mask.shape[1]-8]) # TODO HARDCODED 4
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

        if len(selected_cells) == 0:
            return np.array([])
        
        aligned_fluor = [c.aligned_fluor_mask for c in selected_cells]

        x_size = int(np.median([s.shape[0] for s in aligned_fluor]))
        y_size = int(np.median([s.shape[1] for s in aligned_fluor]))

        resized_cells = self.resize_arr(aligned_fluor, x_size, y_size)
        cell_model = self.create_average(resized_cells)

        return cell_model
    
    def scatter_coords(self, selected_cells):
        
        if not self.spots:
            return np.array([]), np.array([])
        if len(selected_cells) == 0:
            return np.array([]), np.array([])
        
        x_size = np.median([s.aligned_fluor_mask.shape[0] for s in selected_cells])
        y_size = np.median([s.aligned_fluor_mask.shape[1] for s in selected_cells])
        
        x_array = []
        y_array = []
        for cell in selected_cells:
            cy = cell.aligned_fluor_mask.shape[1] / 2 - 0.5
            cx = cell.aligned_fluor_mask.shape[0] / 2 - 0.5
            angle = cell.angle_of_rotation
            for c in cell.spots_coords:
                xc, yc = c
                r_x = cx + math.cos(angle) * (xc - cx) - math.sin(angle) * (yc - cy)
                r_y = cy + math.sin(angle) * (xc - cx) + math.cos(angle) * (yc - cy)
                x_array.append(r_x*(x_size/cell.aligned_fluor_mask.shape[0]))
                y_array.append(r_y*(y_size/cell.aligned_fluor_mask.shape[1]))

        return x_array, y_array

    @staticmethod
    def resize_arr(imgarr:list, xsize:int, ysize:int)->np.ndarray:
        resized = np.zeros((xsize, ysize, len(imgarr)))
        for idx, img in enumerate(imgarr):
            resized[:,:,idx] = resize(img, (xsize, ysize))
        return resized

    @staticmethod
    def create_average(imgarr:np.ndarray)->np.ndarray:
        return np.average(imgarr, axis=2)

    """
    @staticmethod
    def read_xml_spots(xmlfile:str)->Spots:
        # Read xml file into memory
        root = ET.parse(xmlfile).getroot()
        model_child = root[0]  # This is trackmate model object
        allspots = model_child[1]  # This is the spot xml node
        quality = float(root[1][4][0].attrib['value'])

        return Spots(allspots, quality)
    """   
    @staticmethod
    def read_xml_spots(xmlfile:str)->SpotsV2:
        # Read xml file into memory
        root = ET.parse(xmlfile).getroot()
        model_child = root[1]  # This is trackmate model object
        allspots = model_child[1]  # This is the spot xml node

        return SpotsV2(allspots)
    