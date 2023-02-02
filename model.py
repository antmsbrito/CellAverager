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
from cellcycleclassifier import CellCycleClassifier

class Replicate:
    """
    This class controls and stores paths and the model for one replicate of an experimental condition
    """

    def __init__(self,fluor1path, basepath:str, base_type:str, fluor2path=None, membranepath=None)->None:
        
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

        # Check for an xml trackmate file
        if os.path.exists(self.fluor1path.replace('.tif', '.xml')):
            self.xml1path = self.fluor1path.replace('.tif', '.xml')
        else: 
            self.xml1path = None

        if fluor2path:
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
        if self.membpath:
            # Classify cells
            ccc = CellCycleClassifier()
            ccc.classify_cells(imm, cel, "Epifluorescence", False) #TODO Careful with epi! 

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
        """

        cm_object = self.cellmodeler[channel-1]

        if cm_object:
            cells = cm_object.select_cells(minspots, maxspots, cellcycle)
            
            if modeltype == "spot":
                return len(cells), cm_object.build_spot_model(cells)
            elif modeltype == "average":
                return len(cells), cm_object.build_average_model(cells)

        else:
            return 0, None


class ExperimentalCondition:
    """
    This class controls and stores paths and models for all the replicates in an experimental condition
    """

    def __init__(self, root_path, fluor1, fluor2, memb_name, base_name, base_type)->None:

        self.root_path = root_path

        self.fluor1_name = fluor1
        self.fluor2_name = fluor2

        self.memb_name = memb_name

        self.base_name = base_name
        self.base_type = base_type
        
        self.replicates = None

        self.search_root()

    def search_root(self)->None:
        
        replicate_files = list(os.listdir(self.root_path))
        with mp.Pool() as p:
            self.replicates = p.map(self.read_replicate, replicate_files)


    def read_replicate(self, replicate):
        replicate_fullname = os.path.join(self.root_path, replicate)
        if os.path.isdir(replicate_fullname):
            base_path = os.path.join(replicate_fullname, self.base_name + '.tif')
            fluor1_path = os.path.join(replicate_fullname, self.fluor1_name + '.tif')
            fluor2_path = os.path.join(replicate_fullname, self.fluor2_name + '.tif')
            membrane_path = os.path.join(replicate_fullname, self.memb_name + '.tif')

            repli = Replicate(fluor1_path, base_path, self.base_type, fluor2_path, membrane_path)

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
        replicate_models = []

        for repli in self.replicates:
            if repli:
                n, model = repli.buildmodel(channel, modeltype, minspots, maxspots, cellcycle)
                totalN += n
                replicate_models.append(model)

            
        x_size = int(np.median([s.shape[0] for s in replicate_models]))
        y_size = int(np.median([s.shape[1] for s in replicate_models]))

        resized_cells = self.resize_arr(replicate_models, x_size, y_size)
        final_model = np.sum(resized_cells, axis=2) / totalN

        return final_model

    @staticmethod
    def resize_arr(imgarr:list, xsize:int, ysize:int)->np.ndarray:
        resized = np.zeros((xsize, ysize, len(imgarr)))
        for idx, img in enumerate(imgarr):
            resized[:,:,idx] = resize(img, (xsize, ysize))
        return resized