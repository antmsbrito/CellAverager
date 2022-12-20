import os
import base64
from PIL import Image
from io import BytesIO

from matplotlib import pyplot as plt
import numpy as np


def img_to_b64(path2cfp, path2yfp):
    if path2cfp is None:
        cfp_img = np.zeros((30, 30, 3))
    else:
        cfp_img = plt.imread(path2cfp)

    if path2yfp is None:
        yfp_img = np.zeros((30, 30, 3))
    else:
        yfp_img = plt.imread(path2yfp)

    fig = plt.figure()

    ax_cfp = fig.add_subplot(1, 2, 1)
    ax_cfp.set_title("CFP")
    ax_cfp.imshow(cfp_img)
    ax_cfp.axis('off')

    ax_yfp = fig.add_subplot(1, 2, 2)
    ax_yfp.set_title("YFP")
    ax_yfp.imshow(yfp_img)
    ax_yfp.axis('off')

    tmpfile = BytesIO()
    fig.savefig(tmpfile, format='png', bbox_inches='tight', dpi=100)
    encoded = base64.b64encode(tmpfile.getvalue()).decode('utf-8')

    plt.close(fig)

    return encoded


if __name__ == '__main__':

    html = "<!DOCTYPE html><html><body>"
    html += "<h1>Heatmaps 30-11-2022</h1>"

    root = r"C:\Users\António\Documents\30-11-2022_heatmaps"
    for condition in os.listdir(root):

        conditionpath = os.path.join(root, condition)

        if os.path.isdir(conditionpath):

            html += f"<h2>{condition}</h3>"
            yfp = None
            cfp = None
            for img in os.listdir(conditionpath):
                if img.endswith('png') and ('cfp' in img and 'TM_SPOTS' in img):
                    cfp = os.path.join(conditionpath, img)
                elif img.endswith('png') and ('yfp' in img and 'TM_SPOTS' in img):
                    yfp = os.path.join(conditionpath, img)
            html += f'<img src=\'data:image/png;base64,{img_to_b64(cfp, yfp)}\'>'

    html += "</body></html>"
    with open(os.path.join(root, "Heatmap_overview_TM.html"), 'w') as f:
        f.write(html)
