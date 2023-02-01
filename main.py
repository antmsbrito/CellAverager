import os
import sys

from GUI import MainGUI


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
        join_replicates(path, fluornames, (1,))
        join_replicates(path, fluornames, (2,))
        join_replicates(path, fluornames, (3,))
        join_replicates(path, fluornames, (1,2,3))