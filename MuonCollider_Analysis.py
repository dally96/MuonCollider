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

colors  = ["#A4036F","#048ba8","#16db93","#efea5a","#f29e4c","#f96e21","#ff0000","#000000","#ffffff", ]
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
ThetaRes = 9#0.025 #rad
StTimeRes = 30e-3 #ns, st-standard
TimeResScan = [StTimeRes/5,StTimeRes,StTimeRes*2,StTimeRes*5,StTimeRes*10,StTimeRes*50,StTimeRes*100]
tmpNumEntries = tmpPos.size
tmpAbsTimeDiffScan = []
for TimeRes in TimeResScan:
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
    tmpAbsTimeDiffScan.append(tmpAbsTimeDiff)

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

hsTimeDiffScan = []
for TimeRes in TimeResScan:
    hsTime = noBIBTree["sttim"].array()
    hsTime = smear(hsTime,TimeRes)
    hsTimeEx = hsPos/SpeedOfLight
    hsTimeDiff = np.subtract(hsTime,hsTimeEx)
    hsAbsTimeDiff = np.absolute(hsTimeDiff)
    hsTimeDiffScan.append(hsTimeDiff)

###### This whole chunk of code is to unpack all 999 events into 1 single tuple and then find at what index the momentum is greater than 1 GeV.
hsMomCut = []

for event in range(len(hsMom)):
    for hit in hsMom[event]:
        if hit > 1:
            hsMomCut.append(True)
        else:
            hsMomCut.append(False)
            
hsMomCut = np.array(hsMomCut)
hsMomCut = np.where(hsMomCut == True)
hsMomCut = np.array(hsMomCut)
######

hsNumEntries = len(hsMomCut[0])

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

hsAbsTimeDiffScan = []
for hsTimeDiff in hsTimeDiffScan:
    ###### Again this unpacks the 999 tuples into one tuple for our needs
    hsTimeDiffComplete = [y for x in hsTimeDiff for y in x]
    hsThetaDiffComplete = [y for x in hsThetaDiff for y in x]
    ######
    SignalTimeDiff = []
    SignalThetaDiff = []

    for i in hsMomCut[0]:
        SignalTimeDiff.append(hsTimeDiffComplete[i])
        SignalThetaDiff.append(hsThetaDiffComplete[i])
    SignalAbsTimeDiff = np.absolute(SignalTimeDiff)
    SignalAbsThetaDiff = np.absolute(SignalThetaDiff)
    hsAbsTimeDiffScan.append(SignalAbsTimeDiff)
    ######
    

#Create efficiency values for ROC curve

tmpEfficiencies = []
hsEfficiencies  = []

tDeltaValues=np.linspace(0.001,5,10000)
for n, curve in enumerate(hsAbsTimeDiffScan):
    tmpEfficiency = []
    hsEfficiency  = []
    print("Building ROC Curve "+str(n+1)+"/7")
    for v in tDeltaValues: 
        
        tmpDelTimeCut = np.where(tmpAbsTimeDiffScan[n] < v)
        tmpNumKeptEntries = len(tmpDelTimeCut[0])

        hsDelTimeCut = np.where(hsAbsTimeDiffScan[n] < v)
        hsNumKeptEntries = len(hsDelTimeCut[0])

    
        tmpEfficiency.append(1-tmpNumKeptEntries/tmpNumEntries)
        hsEfficiency.append(hsNumKeptEntries/hsNumEntries)
       
    
    tmpEfficiencies.append(tmpEfficiency)
    hsEfficiencies.append(hsEfficiency)
    


fig,ax = plt.subplots()
'''
h=plt.hist2d(tmpThetaDiff, tmpTimeDiff, bins=[30,30])
plt.xlabel("θ-θ$_{ex}$ [rad]"); plt.ylabel("T-ToF [ns]")
plt.scatter(SignalThetaDiff,SignalTimeDiff,s=[[5,]*len(SignalTimeDiff)], c=colors[6])
fig.colorbar(h[3], ax = ax)
# plt.show()

plt.savefig(f"test_{TimeRes}_{ThetaRes}.png")
'''

for n, curve in enumerate(hsEfficiencies):
    hsEfficiency = np.array(hsEfficiencies[n])
    tmpEfficiency = np.array(tmpEfficiencies[n])
    plt.scatter(hsEfficiencies,tmpEfficiencies,c=colors[n])
plt.xlabel("sig Efficiency"); plt.ylabel("1-(bkg Efficiency)")
plt.xlim([0,1]); plt.ylim([0,1])
plt.savefig(f"test_MuonROC_TimeScan_{TimeResScan[0]}-{TimeResScan[6]}.png")
