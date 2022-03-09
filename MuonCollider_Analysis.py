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

from awkward.layout import ListOffsetArray64


def convertNBIBToFrac(x):
    return x/2992.

def smear(arr,sigma):
    #Convert it to a 1D numpy array and perform smearing
    numpy_arr = np.asarray(arr.layout.content)
    smeared_arr = np.random.normal(numpy_arr, sigma)
    #Convert it back to awkward form
    return ak.Array(ListOffsetArray64(arr.layout.offsets, ak.Array(smeared_arr).layout))

SpeedOfLight = constants.c/1e6 # mm/ns

mpl.rcParams['patch.force_edgecolor'] = True
mpl.rcParams['patch.linewidth']       = 0.5
mpl.rcParams['savefig.dpi']           = 300
mpl.rcParams['agg.path.chunksize'] = 10000

colors  = ["#A4036F","#048ba8","#16db93","#efea5a","#f29e4c","#000000","#ffffff"]
markers = ["o","s","D","^","v","<",">","*","X","p"]


# hard scatter only (0 BIB)
noBIBFile = uproot.open(f"lctuple_999.root")
noBIBTree = noBIBFile["LCTupleDefault"]

# BIB only
fullBIBFile = uproot.open("lctuple_2992_trackerSimHits.root")
fullBIBTree =  fullBIBFile["MyLCTuple"]

# w/ full BIB
tmpP = np.array([fullBIBTree["stedp"].array()[0], fullBIBTree["stmox"].array()[0], fullBIBTree["stmoy"].array()[0], fullBIBTree["stmoz"].array()[0]])

tmpX = np.array([fullBIBTree["stedp"].array()[0], fullBIBTree["stpox"].array()[0], fullBIBTree["stpoy"].array()[0], fullBIBTree["stpoz"].array()[0]])

tmpX_Squared = np.square(tmpX[1])
tmpY_Squared = np.square(tmpX[2])
tmpZ_Squared = np.square(tmpX[3])
tmpR = np.sqrt(tmpX_Squared + tmpY_Squared)
tmpPos = np.sqrt(tmpX_Squared + tmpY_Squared + tmpZ_Squared)

tmpPx_Squared = np.square(tmpP[1])
tmpPy_Squared = np.square(tmpP[2])
tmpPz_Squared = np.square(tmpP[3])
tmpPr = np.sqrt(tmpPx_Squared + tmpPy_Squared)
tmpMom = np.sqrt(tmpPx_Squared + tmpPy_Squared + tmpPz_Squared)

##Smearing arrays
ThetaRes = 0 #0.025 #rad
TimeRes = 0 #30e-3 #ns
tmpNumEntries = tmpPos.size
tmpThetaSmr = np.random.normal(0.0,ThetaRes,tmpNumEntries)
tmpTimeSmr =  np.random.normal(0.0,TimeRes, tmpNumEntries)

tmpTheta = np.arctan2(tmpPr,tmpP[3])
tmpTheta = np.add(tmpTheta,tmpThetaSmr)
tmpThetaEx = np.arctan2(tmpR,tmpX[3])
tmpThetaDiff = np.subtract(tmpTheta,tmpThetaEx)
tmpAbsThetaDiff = np.absolute(tmpThetaDiff)

tmpTime = fullBIBTree["sttim"].array()[0]
tmpTime = np.add(tmpTime,tmpTimeSmr)
tmpTimeEx = tmpPos/SpeedOfLight
tmpTimeDiff = np.subtract(tmpTime,tmpTimeEx)
tmpAbsTimeDiff = np.absolute(tmpTimeDiff)

### The next chunks of code are to find at which index the if condition is true.
tmpThetaCut = []

for hit in tmpAbsThetaDiff:
    if hit < 0.1:
        tmpThetaCut.append(True)
    else:
        tmpThetaCut.append(False)

tmpThetaCut = np.array(tmpThetaCut)
tmpThetaCut = np.where(tmpThetaCut == True)
tmpThetaCut = np.array(tmpThetaCut)

tmpTimeCut = []

for hit in tmpAbsTimeDiff:
    if hit < 0.03:
        tmpTimeCut.append(True)
    else:
        tmpTimeCut.append(False)

tmpTimeCut = np.array(tmpTimeCut)
tmpTimeCut = np.where(tmpTimeCut == True)
tmpTimeCut = np.array(tmpTimeCut)


# w/ no BIB, hard scatter only

###The reason we use ak.Array here is because the hard scatter file has 999 events, whereas the BIB file has 1 event, and each of those 999 events have a different number of entries, so it cannot be put into a numpy array. Therefore, each variable extracted from the hard scatter root file is a tuple containing 999 different tuples,each with a different number of events

hsP = ak.Array([noBIBTree["stedp"].array(), noBIBTree["stmox"].array(), noBIBTree["stmoy"].array(), noBIBTree["stmoz"].array()])

hsX = ak.Array([noBIBTree["stedp"].array(), noBIBTree["stpox"].array(), noBIBTree["stpoy"].array(), noBIBTree["stpoz"].array()])


hsX_Squared = np.square(hsX[1])
hsY_Squared = np.square(hsX[2])
hsZ_Squared = np.square(hsX[3])
hsR = np.sqrt(hsX_Squared + hsY_Squared)
hsPos = np.sqrt(hsX_Squared + hsY_Squared + hsZ_Squared)

hsPx_Squared = np.square(hsP[1])
hsPy_Squared = np.square(hsP[2])
hsPz_Squared = np.square(hsP[3])
hsPr = np.sqrt(hsPx_Squared + hsPy_Squared)
hsMom = np.sqrt(hsPx_Squared + hsPy_Squared + hsPz_Squared)

hsTheta = np.arctan2(hsPr,hsP[3])
hsTheta = smear(hsTheta,ThetaRes)
hsThetaEx = np.arctan2(hsR,hsX[3])
hsThetaDiff = np.subtract(hsTheta,hsThetaEx)
hsAbsThetaDiff = np.absolute(hsThetaDiff)

hsTime = noBIBTree["sttim"].array()
hsTime = smear(hsTime,TimeRes)
hsTimeEx = hsPos/SpeedOfLight
hsTimeDiff = np.subtract(hsTime,hsTimeEx)
hsAbsTimeDiff = np.absolute(hsTimeDiff)

###### This whole chunk of code is to unpack all 999 events into 1 single tuple and then find at what index the momentum is greater than 1 GeV.
hsMomCut = []

for event in range(len(hsMom)):
    for hit in hsMom[event]:
        # if hit > 1:
        if hit > 0:
            hsMomCut.append(True)
        else:
            hsMomCut.append(False)

hsMomCut = np.array(hsMomCut)
hsMomCut = np.where(hsMomCut == True)
hsMomCut = np.array(hsMomCut)
######


hsThetaCut = []

for event in range(len(hsAbsThetaDiff)):
    for hit in hsAbsThetaDiff[event]:
        if hit < 0.1:
            hsThetaCut.append(True)
        else:
            hsThetaCut.append(False)

hsThetaCut = np.array(hsThetaCut)
hsThetaCut = np.where(hsThetaCut == True)
hsThetaCut = np.array(hsThetaCut)

hsTimeCut = []

for event in range(len(hsAbsTimeDiff)):
    for hit in range(len(hsAbsTimeDiff[event])):
        if hit < 0.03:
            hsTimeCut.append(True)
        else:
            hsTimeCut.append(False)

hsTimeCut = np.array(hsTimeCut)
hsTimeCut = np.where(hsTimeCut == True)
hsTimeCut = np.array(hsTimeCut)

hsTimeMomCut = np.intersect1d(hsTimeCut, hsMomCut)
hsThetaMomCut = np.intersect1d(hsThetaCut, hsMomCut)
hsCompleteCut = np.intersect1d(hsTimeMomCut, hsThetaMomCut)

###### Again this unpacks the 999 tuples into one tuple for our needs
hsTimeDiffComplete = [y for x in hsTimeDiff for y in x]
hsThetaDiffComplete = [y for x in hsThetaDiff for y in x]
######
SignalTimeDiff = []
SignalThetaDiff = []


for i in hsMomCut[0]:
    SignalTimeDiff.append(hsTimeDiffComplete[i])
    SignalThetaDiff.append(hsThetaDiffComplete[i])

for i in range(20):
    print(SignalTimeDiff[i], hsTimeDiffComplete[i])

print(len(SignalTimeDiff), len(hsTimeDiffComplete) )

######

cmap = plt.cm.viridis.copy()
cmap.set_under(cmap(1e-5),1)

######

def commonStuff(ax,doText=True):
    if doText:
        ax.text(0.05, 0.93, r"$\sqrt{s}=1.5$ TeV Circular Muon Collider",
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes, fontsize=11, c="w")
        ax.text(0.05, 0.88, "MARS15 BIB, CLIC_o3_v14_mod4, Vertex Detector",
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes, fontsize=11, c="w")
        # ax.text(0.05, 0.83, "Vertex Detector",
        #         verticalalignment='bottom', horizontalalignment='left',
        #         transform=ax.transAxes, fontsize=11, c="w")
    plt.tight_layout()
    return


######

fig,ax = plt.subplots()

h=plt.hist2d(tmpThetaDiff, tmpTimeDiff, bins=[200,200], cmap=cmap, rasterized=True)
plt.xlabel(r"$\theta_{\Delta}$ [rad]"); plt.ylabel(r"$t_{\Delta}$ [ns]")
plt.scatter(SignalThetaDiff,SignalTimeDiff, s=2, lw=0.0, ec="k", c="r", alpha=0.02,rasterized=True)
fig.colorbar(h[3], ax = ax, label="BIB Particle Interactions")
commonStuff(ax)

ax.text(0.5,0.05, "Collision Product Interactions", rotation=0,
        verticalalignment='center', horizontalalignment='left',
        multialignment="right", transform=ax.transAxes, fontsize=10, c="w")
plt.plot([0.47],[0.054],"s",c="r",ms=6,mew=0.5,mec="w",alpha=0.8,transform=ax.transAxes)

plt.savefig(f"thetaDelta_tDelta_{TimeRes}_{ThetaRes}.pdf")

######

fig,ax = plt.subplots()

h=plt.hist2d(tmpTime, tmpTimeEx, bins=[200,200], range=[[-0.1, 0.7], [-0.1, 0.7]], cmap=cmap, rasterized=True)
plt.xlabel(r"$t$ [ns]"); plt.ylabel(r"$t_{exp}$ [ns]")
plt.scatter(ak.flatten(hsTime),ak.flatten(hsTimeEx), s=2, lw=0.0, ec="k", c="r", alpha=0.02,rasterized=True)
fig.colorbar(h[3], ax = ax, label="BIB Particle Interactions")
plt.xlim([-0.1,0.7]); plt.ylim([-0.1,0.7])
commonStuff(ax)


ax.text(0.5,0.05, "Collision Product Interactions", rotation=0,
        verticalalignment='center', horizontalalignment='left',
        multialignment="right", transform=ax.transAxes, fontsize=10, c="w")
plt.plot([0.47],[0.054],"s",c="r",ms=6,mew=0.5,mec="w",alpha=1.0,transform=ax.transAxes)

# ax.annotate("Collision Product Interactions",
#             xy=(0.15,0.15), xycoords='data',
#             xytext=(70, -40), textcoords='offset points', fontsize=10, c="r",
#             arrowprops=dict(arrowstyle="simple,head_length=1,head_width=1,tail_width=0.3",
#                             connectionstyle="angle3,angleA=0,angleB=150",
#                             color="r",
#                             ec="w",
#                             ),
#             horizontalalignment='left', verticalalignment='top',
#             )

plt.savefig(f"tEx_t_{TimeRes}_{ThetaRes}.pdf")


######

fig,ax = plt.subplots()

h=plt.hist2d(tmpTheta, tmpThetaEx, bins=[200,200], range=[[0,3.19], [0,3.19]], cmap=cmap, rasterized=True)
plt.xlabel(r"$\theta$ [rad]"); plt.ylabel(r"$\theta_{exp}$ [rad]")
plt.scatter(ak.flatten(hsTheta),ak.flatten(hsThetaEx), s=2, lw=0.0, ec="k", c="r", alpha=0.02,rasterized=True)
fig.colorbar(h[3], ax = ax, label="BIB Particle Interactions")
commonStuff(ax)

ax.text(0.5,0.05, "Collision Product Interactions", rotation=0,
        verticalalignment='center', horizontalalignment='left',
        multialignment="right", transform=ax.transAxes, fontsize=10, c="w")
plt.plot([0.47],[0.054],"s",c="r",ms=6,mew=0.5,mec="w",alpha=0.8,transform=ax.transAxes)

plt.savefig(f"thetaEx_theta_{TimeRes}_{ThetaRes}.pdf")

# plt.show()

######

fig,ax = plt.subplots()

h=plt.hist2d(tmpX[3], tmpTheta, bins=[200,200], cmap=cmap, rasterized=True)
plt.xlabel(r"$z$ [cm]"); plt.ylabel(r"$\theta$ [rad]")
plt.scatter(ak.flatten(hsX[3]),ak.flatten(hsTheta), s=2, lw=0.0, ec="k", c="r", alpha=0.02,rasterized=True)
fig.colorbar(h[3], ax = ax, label="BIB Particle Interactions")
commonStuff(ax,False)

plt.savefig(f"z_theta_{TimeRes}_{ThetaRes}.pdf")

######

fig,ax = plt.subplots()

h=plt.hist2d(tmpX[3], tmpThetaDiff, bins=[200,200], cmap=cmap, rasterized=True)
plt.xlabel(r"$z$ [cm]"); plt.ylabel(r"$\theta_{\Delta}$ [rad]")
plt.scatter(ak.flatten(hsX[3]),ak.flatten(hsThetaDiff), s=2, lw=0.0, ec="k", c="r", alpha=0.02,rasterized=True)
fig.colorbar(h[3], ax = ax, label="BIB Particle Interactions")
commonStuff(ax,False)

plt.savefig(f"z_thetaDelta_{TimeRes}_{ThetaRes}.pdf")

#####

fig,ax = plt.subplots()

h=plt.hist2d(tmpX[3], tmpTime, bins=[200,200], cmap=cmap, rasterized=True)
plt.xlabel(r"$z$ [cm]"); plt.ylabel(r"$t$ [ns]")
plt.scatter(ak.flatten(hsX[3]),ak.flatten(hsTime), s=2, lw=0.0, ec="k", c="r", alpha=0.02,rasterized=True)
fig.colorbar(h[3], ax = ax, label="BIB Particle Interactions")
commonStuff(ax,False)

plt.savefig(f"z_t_{TimeRes}_{ThetaRes}.pdf")


#####

fig,ax = plt.subplots()

h=plt.hist2d(tmpX[3], tmpTimeDiff, bins=[200,200], cmap=cmap, rasterized=True)
plt.xlabel(r"$z$ [cm]"); plt.ylabel(r"$t_{\Delta}$ [ns]")
# plt.scatter(ak.flatten(hsX[3]),ak.flatten(hsTimeDiff), s=2, lw=0.5, ec="k", c=colors[6], alpha=0.02,rasterized=True)
plt.scatter(ak.flatten(hsX[3]),ak.flatten(hsTimeDiff), s=2, lw=0.0, ec="k", c="r", alpha=0.02,rasterized=True)
fig.colorbar(h[3], ax = ax, label="BIB Particle Interactions")
commonStuff(ax,False)

plt.savefig(f"z_tDelta_{TimeRes}_{ThetaRes}.pdf")