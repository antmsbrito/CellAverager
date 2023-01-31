"""
Classes to help read trackmate xml's
Only deals with spot detection results
Antonio Brito @ ITQB BCB 2022
"""

import numpy as np


class Spots:

    def __init__(self, all_spots_node, minquality):

        self.n_spots = int(all_spots_node.attrib['nspots'])
        self.minquality = minquality

        self.allspots = []
        self.quality = []

        for frame in all_spots_node:
            for spot in frame:
                self.allspots.append(spot)
                self.quality.append(float(spot.attrib['QUALITY']))
        self.quality = np.array(self.quality)

        # Sanity check
        assert self.n_spots == len(self.allspots)

    def get_position(self, sp):

        # x axis is horizontal, right to left
        x = float(sp.attrib['POSITION_X'])  # px size!? TODO currently this has to be manually changed case by case
        # y axis is vertical UP to DOWN
        y = float(sp.attrib['POSITION_Y'])  # px size!? TODO currently this has to be manually changed case by case
        return x, y

    def filterbox(self, box, align):
        # THESE ARE INDICES
        # x is lines => Y axis
        # y is columns => X axis
        x0, y0, x1, y1 = box

        count = 0
        coords = []
        for sp in self.allspots:

            if float(sp.attrib['QUALITY']) < self.minquality:
                continue

            xc, yc = self.get_position(sp)
            xc, yc = xc + (align[1]), yc + (align[0])

            if x0 < yc < x1 and y0 < xc < y1:
                # the coords of the spot in the local reference frame (y0,x0) = (0,0)
                coords.append((y1 - xc, x1 - yc))
                count += 1

        return count, coords
