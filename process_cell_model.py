import numpy as np
import matplotlib as mpl
from skimage.io import imread, imsave
from skimage.filters import threshold_isodata
from skimage.exposure import rescale_intensity
from skimage.color import gray2rgb
from skimage.util import img_as_uint
from matplotlib import pyplot as plt
from tkinter import filedialog as fd


def assign_color(filtered, color, cmap):
    # todo change color array input to shape tuple
    norm = mpl.colors.Normalize(vmin=np.min(filtered[np.nonzero(filtered)]), vmax=np.max(filtered))
    for i in range(filtered.shape[0]):
        for ii in range(filtered.shape[1]):
            px = filtered[i, ii]
            if px == 0:
                color[i, ii] = cmap(0)[:-1]
            else:
                rgba = cmap(norm(px))
                color[i, ii] = (rgba[0], rgba[1], rgba[2])
    return color


def model2color(path2model, savepath):
    model = imread(path2model)

    mask = model > threshold_isodata(model)

    filtered = model * mask

    colormap = mpl.cm.get_cmap("coolwarm")
    color_img = np.zeros(np.shape(gray2rgb(filtered)))

    color_img = assign_color(filtered, color_img, colormap)

    imsave(savepath, color_img)
