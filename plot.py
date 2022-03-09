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
(0,0),
(0.02,0.01),
(0.02,0.05),
(0.05,0.01),
(0.05,0.05),
]

starttimecut=0
# endtimecut=starttimecut+200
endtimecut=-1

fig, ax = plt.subplots(nrows=1, ncols=1)

# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/zoom_inset_axes.html

for imap,(TimeRes,ThetaRes) in enumerate(pairs):

    with open(f"fromLukeTest/EffDict_{TimeRes}ns_{ThetaRes}rad.pickle","rb") as f:
    # with open(sys.argv[1],"rb") as f:
        object_file = pickle.load(f)

    [sigEffsDict,bgEffsDict] = object_file

    print(bgEffsDict.keys())
    print(sigEffsDict.keys())

    colors = cmaps[imap](np.linspace(0.0, 1.0, len(bgEffsDict)-starttimecut ) ) 
    icolor=0
    timecutlist = list(sigEffsDict.keys() )
    start,end = 0,-1
    for i,timecut in enumerate(timecutlist[starttimecut:endtimecut]):
        sigEffsDict[timecutlist[i]] = sigEffsDict[timecutlist[i]][start:end]
        bgEffsDict[timecutlist[i]] = bgEffsDict[timecutlist[i]][start:end]
        plt.fill(
            np.append(
                sigEffsDict[timecutlist[i]],
                sigEffsDict[timecutlist[i+1]][::-1]
                ),
            np.append(
                bgEffsDict[timecutlist[i]],
                bgEffsDict[timecutlist[i+1]][::-1]
                ),
            c=colors[icolor],alpha=1.0
          )
        icolor+=1

ax.set_xlabel("Collision Product Efficiency")
ax.set_ylabel("BIB Efficiency")
ax.set_yscale('log')    
# ax.set_xscale('log')    






# inset axes....
axins = ax.inset_axes([0.1, 0.65, 0.35, 0.33])

for imap,(TimeRes,ThetaRes) in enumerate(pairs):

    with open(f"fromLukeTest/EffDict_{TimeRes}ns_{ThetaRes}rad.pickle","rb") as f:
        object_file = pickle.load(f)

    [sigEffsDict,bgEffsDict] = object_file
    colors = cmaps[imap](np.linspace(0.0, 0.8, len(bgEffsDict)-starttimecut ) ) 
    icolor=0
    timecutlist = list(sigEffsDict.keys() )
    for i,timecut in enumerate(timecutlist[starttimecut:-1]):
        axins.fill(
            np.append(sigEffsDict[timecutlist[i]],sigEffsDict[timecutlist[i+1] ][::-1]),
            np.append(bgEffsDict[timecutlist[i]],bgEffsDict[timecutlist[i+1]][::-1]),
            c=colors[icolor],alpha=1.0
          )
        icolor+=1

# sub region of the original image
x1, x2, y1, y2 = 0.96,1.005,8e-3,1.2e-1
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_yscale('log')    

rect,c = ax.indicate_inset_zoom(axins, edgecolor="black",alpha=0.8,lw=0.5)
plt.setp(c, lw=0.3,alpha=0.4)


imap=2
with open(f"fromLukeTest/EffDict_0.02ns_0.05rad.pickle","rb") as f:
    object_file = pickle.load(f)

[sigEffsDict,bgEffsDict] = object_file
colors = cmaps[imap](np.linspace(0.0, 0.8, len(bgEffsDict)-starttimecut ) ) 
icolor=0
timecutlist = list(sigEffsDict.keys() )
icut = 10
jcut = 6
x,y=sigEffsDict[timecutlist[icut]][jcut], bgEffsDict[timecutlist[icut]][jcut]
thetaCut = np.linspace(0.5*ThetaRes,5*ThetaRes,20)[jcut]
axins.plot( x,y, "ow",mew=0.5,mec="k" )
axins.annotate(
    rf'$|t_\Delta|<{timecutlist[icut]:0.3f}$, $|\theta_\Delta|<{thetaCut:0.3f}$',
    xy=(x,y),  xycoords='data',
            xytext=(1.4, 0.35), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="simple,head_width=0.9,head_length=0.9",
                            connectionstyle="arc3,rad=0.3",facecolor="k"),
            horizontalalignment='left', verticalalignment='top',
            )


# legend
def drawLegendEntry(x,y,label,cmap):
    ax.text(x,y, label,
            verticalalignment='center', horizontalalignment='left',
            multialignment="left", transform=ax.transAxes, fontsize=10, c="k")
    plt.plot([x-0.03],[y+0.004],"s",c=cmap(0.2),ms=6,mew=0.5,mec="k",alpha=0.8,transform=ax.transAxes)

drawLegendEntry(0.65,0.18+0.05,r"$\sigma(t)=0$ ps, $\sigma(\theta)=0$ rad",cmaps[0])
drawLegendEntry(0.65,0.14+0.05,r"$20$ ps, $0.01$ rad",cmaps[1])
drawLegendEntry(0.65,0.10+0.05,r"$20$ ps, $0.05$ rad",cmaps[2])
drawLegendEntry(0.65,0.06+0.05,r"$50$ ps, $0.01$ rad",cmaps[3])
drawLegendEntry(0.65,0.02+0.05,r"$50$ ps, $0.05$ rad",cmaps[4])








plt.tight_layout()

plt.savefig(f"ROC.pdf")
plt.show()







