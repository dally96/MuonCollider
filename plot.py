import uproot
import random
import awkward as ak
#from uproot_methods import TLorentzVectorArray
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys, os
from scipy import constants
import math
from math import pi

import pickle
import sys
from awkward.layout import ListOffsetArray64

# TimeResList = [0.03,0.05,0.1]
# ThetaResList = [0.015,0.025,0.035]

import matplotlib.cm as cm

cmaps = [
    # cm.get_cmap("spring"),
    # cm.get_cmap("summer"),
    # cm.get_cmap("autumn"),
    # cm.get_cmap("winter"),


    cm.get_cmap("Purples_r"),
    cm.get_cmap("Blues_r"),
    cm.get_cmap("Greens_r"),
    cm.get_cmap("Oranges_r"),
    cm.get_cmap("Reds_r"),
    cm.get_cmap("Greys_r"),
    ]

pairs = [
(0.1,0.035),
(0.05,0.025),
(0.03,0.015),
# (0.1,0.035),
# (0.1,0.015),
]


fig, ax = plt.subplots(nrows=1, ncols=1)


for imap,(TimeRes,ThetaRes) in enumerate(pairs):

    with open(f"effs_{TimeRes}_{ThetaRes}.pkl","rb") as f:

    # with open(sys.argv[1],"rb") as f:
        object_file = pickle.load(f)


    skip=15
    [sigEffsDict,bgEffsDict] = object_file

    colors = cmaps[imap](np.linspace(0.0, 0.8, len(bgEffsDict)-skip ) ) 
    icolor=0
    for i in range(len(bgEffsDict)-1)[skip:]:
        plt.fill(
            np.append(sigEffsDict[i],sigEffsDict[i+1][::-1]),
            np.append(bgEffsDict[i],bgEffsDict[i+1][::-1]), c=colors[icolor]
          )
        icolor+=1

plt.savefig(f"ROC.pdf")

plt.show()







