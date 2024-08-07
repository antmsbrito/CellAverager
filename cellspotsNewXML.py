"""
Classes to help read trackmate xml's
Only deals with spot detection results
Antonio Brito @ ITQB BCB 2022
"""

import numpy as np


class Spots:

    def __init__(self, all_spots_node):

        self.n_spots = int(all_spots_node.attrib['nspots'])

        self.allspots = []

        for frame in all_spots_node:
            for spot in frame:
                self.allspots.append(spot)
        # Sanity check
        assert self.n_spots == len(self.allspots)

    def get_position(self, sp):

        # x axis is horizontal, right to left
        x = float(sp.attrib['POSITION_X']) / 0.08  # px size!? TODO currently this has to be manually changed case by case
        # y axis is vertical UP to DOWN
        y = float(sp.attrib['POSITION_Y']) / 0.08  # px size!? TODO currently this has to be manually changed case by case

        return x, y

    def filterbox(self, box, align):
        # THESE ARE INDICES
        # x is lines => Y axis
        # y is columns => X axis
        x0, y0, x1, y1 = box
        
        x0 = x0 + 4
        y0 = y0 + 4
        
        x1 = x1 - 4
        y1 = y1 - 4
        
        count = 0
        coords = []
        for sp in self.allspots:

            xc, yc = self.get_position(sp)
            xc, yc = xc + (align[1]), yc + (align[0])

            if x0 < yc < x1 and y0 < xc < y1:
                # the coords of the spot in the local reference frame (y0,x0) = (0,0)
                coords.append((y1 - xc, x1 - yc))
                count += 1

        return count, coords

    def filterlabel(self,box,labelid,labelimg):
        
        x0, y0, x1, y1 = box
        # TODO HARDCODED 4
        x0 = x0 + 4
        y0 = y0 + 4
        
        x1 = x1 - 4
        y1 = y1 - 4
        
        count = 0
        coords = []
        for sp in self.allspots:
            
            xc,yc = self.get_position(sp)
            ixc = round(xc)
            iyc = round(yc)
            #todo align
            if x0 < yc < x1 and y0 < xc < y1:
                if labelimg[iyc,ixc]==labelid:
                    coords.append((y1 - xc, x1 - yc))
                    count += 1
                    
        return count,coords