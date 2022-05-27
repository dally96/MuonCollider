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


    cm.get_cmap("Greys_r"),
    cm.get_cmap("Reds_r"),
    cm.get_cmap("Greens_r"),
    cm.get_cmap("Blues_r"),
    cm.get_cmap("Purples_r"),
    cm.get_cmap("Greys_r"),
    ]

pairs = [
(0.0,0.0),

(0.02,0.1),
(0.05,0.1),
(0.1,0.1),
(0.05,2),
]

# starttimecut=0
# # endtimecut=starttimecut+200
# endtimecut=-1

fig, ax = plt.subplots(nrows=1, ncols=1)

# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/zoom_inset_axes.html

for imap,(TimeRes,ThetaRes) in enumerate(pairs):

    # with open(f"fromLarry_052722/EffDict_{TimeRes}_{ThetaRes}.pickle","rb") as f:
    with open(f"EffDict_{TimeRes}_{ThetaRes}.pickle","rb") as f:
    # with open(sys.argv[1],"rb") as f:
        object_file = pickle.load(f)

    [sigEffsDict,bgEffsDict] = object_file

    # print(bgEffsDict.keys())
    # print(sigEffsDict.keys())

    colors = cmaps[imap](np.linspace(0.15, 0.8, len(bgEffsDict) ) )
    icolor=0
    timecutlist = list(sigEffsDict.keys() )
    start,end = 0,-1
    if (TimeRes,ThetaRes)==(0,0):
        timecutlist = timecutlist[1:]
        # colors = cmaps[imap](np.linspace(-2.0, 0.8, len(bgEffsDict) ) )

    for i,timecut in enumerate(timecutlist[:-1]):
        sigEffsDict[timecutlist[i]] = sigEffsDict[timecutlist[i]][start:end]
        bgEffsDict[timecutlist[i]] = bgEffsDict[timecutlist[i]][start:end]
        plt.fill(
            np.append(
                sigEffsDict[timecutlist[i]][:],
                sigEffsDict[timecutlist[i+1]][-2::-1]
                ),
            np.append(
                bgEffsDict[timecutlist[i]][:],
                bgEffsDict[timecutlist[i+1]][-2::-1]
                ),
            c=colors[icolor],alpha=1.0
          )
        icolor+=1

ax.set_xlabel("Collision Product Efficiency")
ax.set_ylabel("BIB Efficiency")
ax.set_yscale('log')

# ax.set_xlim(0.7, 1.01)
# ax.set_ylim(5e-4, 1.5e0)
# ax.set_xscale('log')



# ax.set_xlim(0.5, 1.02)
# ax.set_ylim(5e-3, 1.5e0)
ax.set_xlim(0.0, 1.02)
ax.set_ylim(10e-4, 5e0)



# inset axes....
axins = ax.inset_axes([0.1, 0.65, 0.35, 0.33])

for imap,(TimeRes,ThetaRes) in enumerate(pairs):

    # with open(f"fromLarry_052722/EffDict_{TimeRes}_{ThetaRes}.pickle","rb") as f:
    with open(f"EffDict_{TimeRes}_{ThetaRes}.pickle","rb") as f:
        object_file = pickle.load(f)

    [sigEffsDict,bgEffsDict] = object_file

    colors = cmaps[imap](np.linspace(0.15, 0.8, len(bgEffsDict) ) )
    icolor=0
    timecutlist = list(sigEffsDict.keys() )
    start,end = 0,-1
    if (TimeRes,ThetaRes)==(0,0):
        timecutlist = timecutlist[1:]
        # colors = cmaps[imap](np.linspace(-2.0, 0.8, len(bgEffsDict) ) )

    for i,timecut in enumerate(timecutlist[:-1]):
        sigEffsDict[timecutlist[i]] = sigEffsDict[timecutlist[i]][start:end]
        bgEffsDict[timecutlist[i]] = bgEffsDict[timecutlist[i]][start:end]
        axins.fill(
            np.append(
                sigEffsDict[timecutlist[i]][:],
                sigEffsDict[timecutlist[i+1]][-2::-1]
                ),
            np.append(
                bgEffsDict[timecutlist[i]][:],
                bgEffsDict[timecutlist[i+1]][-2::-1]
                ),
            c=colors[icolor],alpha=1.0
          )
        icolor+=1

# sub region of the original image
x1, x2, y1, y2 = 0.98,1.00,9e-2,1.1
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_yscale('log')

axins.plot([0.99,0.99], [1e-3,1e-0],    "--",c="k",lw=1.5)
# axins.plot([0.98,0.99], [1.8e-2,1.8e-2],"--",c="k",lw=1.5)
axins.plot([0.98,0.99], [y:=0.14,y],"--",c="k",lw=1.5)
axins.plot([0.98,0.99], [y:=0.55,y],"--",c="k",lw=1.5)


rect,c = ax.indicate_inset_zoom(axins, edgecolor="black",alpha=0.8,lw=0.5)
plt.setp(c, lw=0.3,alpha=0.4)






# imap=1
# with open(f"fromLarry/EffDict_0.02_0.01.pickle","rb") as f:
#     object_file = pickle.load(f)
# [sigEffsDict,bgEffsDict] = object_file
# colors = cmaps[imap](np.linspace(0.15, 0.8, len(bgEffsDict) ) )
# icolor=0
# timecutlist = list(sigEffsDict.keys() )
# icut = 30
# jcut = 8
# x,y=sigEffsDict[timecutlist[icut]][jcut], bgEffsDict[timecutlist[icut]][jcut]
# thetaCut = np.linspace(0.5*ThetaRes,10*ThetaRes,20)[jcut] #fix
# # axins.plot( x,y, "ow",mew=0.5,mec="k" )
# axins.annotate(
#     rf'$|t_\Delta|<{timecutlist[icut]:0.3f}$, $|\theta_\Delta|<{thetaCut:0.3f}$',
#     xy=(x,y),  xycoords='data',
#             xytext=(1.1, 0.25), textcoords='axes fraction',
#             arrowprops=dict(arrowstyle="simple,head_width=0.5,head_length=0.5,tail_width=0.1",
#                             connectionstyle="arc3,rad=-0.1",facecolor="k"),
#             horizontalalignment='left', verticalalignment='top',
#             )


# imap=4
# with open(f"fromLarry/EffDict_0.05_0.05.pickle","rb") as f:
#     object_file = pickle.load(f)
# [sigEffsDict,bgEffsDict] = object_file
# colors = cmaps[imap](np.linspace(0.15, 0.8, len(bgEffsDict) ) )
# icolor=0
# timecutlist = list(sigEffsDict.keys() )
# icut = 63
# jcut = 19
# x,y=sigEffsDict[timecutlist[icut]][jcut], bgEffsDict[timecutlist[icut]][jcut]
# print(x,y)
# thetaCut = np.linspace(0.5*ThetaRes,10*ThetaRes,20)[jcut]
# # axins.plot( x,y, "ow",mew=0.5,mec="k" )
# axins.annotate(
#     rf'$|t_\Delta|<{timecutlist[icut]:0.3f}$, $|\theta_\Delta|<{thetaCut:0.3f}$',
#     xy=(x,y),  xycoords='data',
#             xytext=(1.2, 0.65), textcoords='axes fraction',
#             arrowprops=dict(arrowstyle="simple,head_width=0.5,head_length=0.5,tail_width=0.1",
#                             connectionstyle="arc3,rad=-0.05",facecolor="k"),
#             horizontalalignment='left', verticalalignment='top',
#             )










# legend
def drawLegendEntry(x,y,label,cmap):
    ax.text(x,y, label,
            verticalalignment='center', horizontalalignment='left',
            multialignment="left", transform=ax.transAxes, fontsize=10, c="k")
    plt.plot([x-0.03],[y+0.004],"s",c=cmap(0.2),ms=6,mew=0.5,mec="k",alpha=1.0,transform=ax.transAxes)



# (0.02,0.1),
# (0.05,0.1),
# (0.1,0.1),
# (0.05,2),

drawLegendEntry(0.6,0.02+0.03,r"$0$ ps, $0$ rad",cmaps[0])
drawLegendEntry(0.6,0.06+0.03,r"$20$ ps, $0.1$ rad",cmaps[1])
drawLegendEntry(0.6,0.10+0.03,r"$50$ ps, $0.1$ rad",cmaps[2])
drawLegendEntry(0.6,0.14+0.03,r"$100$ ps, $0.1$ rad",cmaps[3])
drawLegendEntry(0.6,0.18+0.03,r"$\sigma(t)=50$ ps, $\sigma(\theta)=2$ rad",cmaps[4])



# (0.0,0.0),
# (0.01,0.01),
# (0.02,0.01),
# (0.05,0.01),
# # (0.02,0.05),
# (0.05,0.05),

# drawLegendEntry(0.6,0.18+0.03,r"$\sigma(t)=0$ ps, $\sigma(\theta)=0$ rad",cmaps[0])
# drawLegendEntry(0.6,0.14+0.03,r"$20$ ps, $0.01$ rad",cmaps[1])
# drawLegendEntry(0.6,0.10+0.03,r"$20$ ps, $0.05$ rad",cmaps[2])
# drawLegendEntry(0.6,0.06+0.03,r"$50$ ps, $0.01$ rad",cmaps[3])
# drawLegendEntry(0.6,0.02+0.03,r"$50$ ps, $0.05$ rad",cmaps[4])


# ax.text(0.98, 0.94, r"$\sqrt{s}=1.5$ TeV Circular Muon Collider",
#         verticalalignment='bottom', horizontalalignment='right',
#         transform=ax.transAxes, fontsize=10, c="k")
ax.text(0.98, 0.84, r"$\sqrt{s}=1.5$ "+"TeV Circular Muon Collider\nMARS15 BIB, CLIC_o3_v14_mod4\nVertex Detector",
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes, fontsize=10, c="k",backgroundcolor='1')





plt.tight_layout()

plt.savefig(f"ROC.pdf")
plt.show()







