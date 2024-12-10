import os
import multiprocessing as mp

import numpy as np
from skimage.transform import resize

from cellmodeler import CellModeler
from cellaligner import CellAligner

from images import ImageManager
from parameters import MaskParameters, RegionParameters, CellParameters, ParametersManager
from segments import SegmentsManager
from cells import CellManager

from cellcycleclassifier_v2 import CellCycleClassifier
#from cellcycleclassifier import CellCycleClassifier

import pandas as pd

class Replicate:
    """
    This class controls and stores paths and the model for one replicate of an experimental condition
    """

    def __init__(self,fluor1path, basepath:str, base_type:str, fluor2path=None, membranepath=None, dnapath=None)->None:
        
        self.fluor1path = fluor1path
        self.basepath = basepath
        self.base_type = base_type
        
        if os.path.exists(fluor2path):
            self.fluor2path = fluor2path
        else:
            self.fluor2path = None

        if os.path.exists(membranepath):
            self.membpath = membranepath
        else:
            self.membpath = None
            
        if os.path.exists(dnapath):
            self.dnapath = dnapath
        else:
            self.dnapath = None

        # Check for an xml trackmate file
        if os.path.exists(self.fluor1path.replace('.tif', '.xml')):
            self.xml1path = self.fluor1path.replace('.tif', '.xml')
        else: 
            self.xml1path = None

        if self.fluor2path:
            if os.path.exists(self.fluor2path.replace('.tif', '.xml')):
                self.xml2path = self.fluor2path.replace('.tif', '.xml')
            else: 
                self.xml2path = None

        # Store ehooke base classes 
        self.imagemanager = [None, None]
        self.cellmanager = [None, None]
        self.cellmodel = [None, None]

        # Store CellAverager classes
        self.cellaligner = [None, None]
        self.cellmodeler = [None, None]
        
        self.run_eHooke(channel=1)
        if fluor2path:
            self.run_eHooke(channel=2)
        
        self.build_CA()
        
    def run_eHooke(self, channel:int)->None:
        """
        Runs eHooke ONCE for the desired fluorescence channel
        """

        if channel==1 or not self.fluor2path:
            fluorescence = self.fluor1path
        elif channel==2 and self.fluor2path:
            fluorescence = self.fluor2path

        par = ParametersManager()

        # Majority of default parameters work. Change some based on the type of base img used
        if self.base_type == "Phase":
            par.imageloaderparams.mask_fill_holes = True
        elif self.base_type == "BF":
            par.imageloaderparams.mask_algorithm = 'StarDist_BF'
        elif self.base_type == "Membrane":
            par.imageloaderparams.mask_algorithm = 'StarDist'
            
        par.imageloaderparams.auto_align = False
        par.cellprocessingparams.cell_force_merge_below = 0
        par.imageloaderparams.mask_dilation = 1

        # Load base img and compute mask
        imm = ImageManager()
        imm.load_base_image(self.basepath, par.imageloaderparams)
        imm.compute_base_mask(par.imageloaderparams)
        imm.compute_mask(par.imageloaderparams)

        # Check if membrane image exists
        # If there is membrane, put membrane as main fluoresence
        if self.membpath:
            imm.load_fluor_image(self.membpath, par.imageloaderparams)
            imm.load_option_image(fluorescence, par.imageloaderparams)
        else:
            imm.load_fluor_image(fluorescence, par.imageloaderparams)

        seg = SegmentsManager()
        seg.compute_segments(par.imageprocessingparams, imm)

        # QC
        seg.save_labels(filename=os.path.join(os.path.dirname(self.basepath), 'Labels'))

        cel = CellManager(par)
        cel.compute_cells(par, imm, seg);
        cel.process_cells(par.cellprocessingparams, imm);
        cel.overlay_cells(imm);

        # If there is membrane, also do cell cycle
        if self.dnapath:
            # Classify cells
            ccc = CellCycleClassifier(imm.fluor_image,self.dnapath,imm,cel)
            ccc.classify_cells("Epifluorescence") #TODO Careful with epi! 

        # Finally, store the cell manager object

        self.cellmanager[channel-1] = cel
        self.imagemanager[channel-1] = imm

    def build_CA(self,)->None:
        """
        Builds cell model objects for both fluorescences (if possible)
        """

        if self.membpath:
            self.cellaligner[0] = CellAligner(self.cellmanager[0], "Optional")
            self.cellmodeler[0] = CellModeler(self.cellmanager[0], self.imagemanager[0], self.xml1path)
            if self.fluor2path:
                self.cellaligner[1] = CellAligner(self.cellmanager[1], "Optional")
                self.cellmodeler[1] = CellModeler(self.cellmanager[1], self.imagemanager[1], self.xml2path)
        else:
            self.cellaligner[0] = CellAligner(self.cellmanager[0], "Main")
            self.cellmodeler[0] = CellModeler(self.cellmanager[0], self.imagemanager[0], self.xml1path)
            if self.fluor2path:
                self.cellaligner[1] = CellAligner(self.cellmanager[1], "Main")
                self.cellmodeler[1] = CellModeler(self.cellmanager[1], self.imagemanager[1], self.xml2path)    
        
    def buildmodel(self, channel=1, modeltype='spot', minspots=1, maxspots=np.inf, cellcycle=(0,1,2,3)):
        """
        Builds a model of a channel with selected cells
        Returns number of cells used and the model in question
        selected_cells, total_cells, n_spots, model
        """

        cm_object = self.cellmodeler[channel-1]

        if cm_object:
            total_cells = cm_object.total_cells
            cells, n_spots = cm_object.select_cells(minspots, maxspots, cellcycle)
            selected_cells = len(cells)
            
            if modeltype == "spot":
                return selected_cells, total_cells, n_spots, cm_object.build_spot_model(cells)
            elif modeltype == "average":
                return selected_cells, total_cells, 0, cm_object.build_average_model(cells)

        else:
            return 0, 0, 0, None
    
    def get_spot_coords(self, channel=1, minspots=1, maxspots=np.inf, cellcycle=(0,1,2,3)):
        
        cm_object = self.cellmodeler[channel-1]

        if cm_object:
            cells, n_spots = cm_object.select_cells(minspots, maxspots, cellcycle)
            return cm_object.scatter_coords(cells)

        else:
            return None,None
        
class ExperimentalCondition:
    """
    This class controls and stores paths and models for all the replicates in an experimental condition
    """

    def __init__(self, root_path, fluor1, fluor2, memb_name, base_name, base_type, dna_name)->None:

        self.root_path = root_path

        self.fluor1_name = fluor1
        self.fluor2_name = fluor2

        self.memb_name = memb_name
        
        self.dna_name = dna_name

        self.base_name = base_name
        self.base_type = base_type
        
        self.replicates = None

        self.search_root()

    def search_root(self)->None:
        
        replicate_files = list(os.listdir(self.root_path))
        #replicate_files = list(filter(str.isdigit,replicate_files)) # TODO BE CAREFUL
        #self.replicates = map(self.read_replicate, replicate_files)
        with mp.Pool() as p:
            self.replicates = p.map(self.read_replicate, replicate_files)
        

    def read_replicate(self, replicate):
        replicate_fullname = os.path.join(self.root_path, replicate)
        if os.path.isdir(replicate_fullname):
            base_path = os.path.join(replicate_fullname, self.base_name + '.tif')
            fluor1_path = os.path.join(replicate_fullname, self.fluor1_name + '.tif')
            fluor2_path = os.path.join(replicate_fullname, self.fluor2_name + '.tif')
            membrane_path = os.path.join(replicate_fullname, self.memb_name + '.tif')
            dna_path = os.path.join(replicate_fullname, self.dna_name + '.tif')

        
            repli = Replicate(fluor1_path, base_path, self.base_type, fluor2_path, membrane_path, dna_path)
            return repli
        else:
            return None
        

    def askformodel(self, channel=1, modeltype='spot', minspots=1, maxspots=np.inf, cellcycle=(0,1,2,3)):

        if channel not in [1,2]:
            raise ValueError(f"Only 1 or 2 are valid channels. You provided {channel}")

        if modeltype not in ['spot', 'average']:
            raise ValueError(f"Only 'spot' or 'average' are valid models. You provided {modeltype}")

        if minspots>maxspots:
            minspots, maxspots = maxspots, minspots
            

        totalN = 0
        totalCells = 0
        totalSpots = 0
        replicate_models = []
        
        for repli in self.replicates:
            if repli:
                try:
                    selected_cells, total_cells, n_spots, model = repli.buildmodel(channel, modeltype, minspots, maxspots, cellcycle)
                    if selected_cells > 0:
                        totalN += selected_cells
                        totalCells += total_cells
                        totalSpots += n_spots
                        replicate_models.append(model*selected_cells)
                except Exception as e:
                    print(e, repli.fluor1path)
        
        x_size = int(np.median([s.shape[0] for s in replicate_models]))
        y_size = int(np.median([s.shape[1] for s in replicate_models]))

        resized_cells = self.resize_arr(replicate_models, x_size, y_size)
        final_model = np.sum(resized_cells, axis=2) / totalN

        return totalN, totalCells, totalSpots, final_model
    
    def askforscatter(self,channel=1, minspots=1, maxspots=np.inf, cellcycle=(0,1,2,3)):
        
        x = np.array([])
        y = np.array([])
        for repli in self.replicates:
            if repli:
                xcoords,ycoords = repli.get_spot_coords(channel=channel, minspots=minspots, maxspots=maxspots, cellcycle=cellcycle)
                x = np.concatenate([x,xcoords])
                y = np.concatenate([y,ycoords])
        return x,y
    
    def askforscatterpercell(self,channel=1,cellcycle=(0,1,2,3)):
        
        coords_per_cell = {}
        for idxrepli, repli in enumerate(self.replicates):
            if repli:
                for key,cell in repli.cellmanager[channel-1].cells.items():
                    if cell.stats["Cell Cycle Phase"] in cellcycle:
                        coords_per_cell[str(key)+f"_{idxrepli}"] = cell.spots_coords #TODO these are the coords in the local reference frame of the cell bounding box without the 4px border
        
        return coords_per_cell 
    
    def askfortable(self,channel=1):
        
        fov = []
        cellid = []
        nspots = []
        ccp = []
        spot_coords = []
        
        for idx,repli in enumerate(self.replicates):
            if repli:
                cm = repli.cellmanager[channel-1]
                for cellkey in cm.cells:
                    fov.append(os.path.basename(os.path.dirname(repli.fluor1path)))
                    cellid.append(cellkey)
                    nspots.append(cm.cells[cellkey].spots)
                    ccp.append(cm.cells[cellkey].stats['Cell Cycle Phase'])
                    spot_coords.append(cm.cells[cellkey].spots_coords)
                    
        
        return pd.DataFrame(data={'FoV':fov,'CellID':cellid,f'Nspots_{channel}':nspots,'CCP':ccp, f'SpotCoords_{channel}':spot_coords})
        
    @staticmethod
    def resize_arr(imgarr:list, xsize:int, ysize:int)->np.ndarray:
        resized = np.zeros((xsize, ysize, len(imgarr)))
        for idx, img in enumerate(imgarr):
            resized[:,:,idx] = resize(img, (xsize, ysize))
        return resized