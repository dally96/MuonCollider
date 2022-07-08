import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle

with open(f"DiHiggs_ZREThExDelTDel_ECut.pickle","rb") as f:
    DiHiggs = pickle.load(f)

with open(f"mmTOhmm_ZREThExDelTDel_ECut.pickle","rb") as f:
    mmTOhmm = pickle.load(f)

with open(f"mmTOhvv_ZREThExDelTDel_ECut.pickle","rb") as f:
    mmTOhvv = pickle.load(f)

with open(f"mmTOh_ZREThExDelTDel_ECut.pickle","rb") as f:
    mmTOh = pickle.load(f)

with open(f"BIB_ZREThExDelTDel_ECut.pickle","rb") as f:
    BIB = pickle.load(f)

with open(f"DiHiggs_EffDicts.pickle","rb") as f:
    DiHiggs_Eff = pickle.load(f)

with open(f"mmTOhmm_EffDicts.pickle","rb") as f:
    mmTOhmm_Eff = pickle.load(f)

with open(f"mmTOhvv_EffDicts.pickle","rb") as f:
    mmTOhvv_Eff = pickle.load(f)

with open(f"mmTOh_EffDicts.pickle","rb") as f:
    mmTOh_Eff = pickle.load(f)

ThMPlots = []

ThPPlots = []
ThPPlots.append(DiHiggs[4])
ThPPlots.append(mmTOhmm[4])
ThPPlots.append(mmTOhvv[4])
ThPPlots.append(mmTOh[4])
ThPPlots.append(BIB[4])

DelThPlots = []
DelThPlots.append(DiHiggs[5])
DelThPlots.append(mmTOhmm[5])
DelThPlots.append(mmTOhvv[5])
DelThPlots.append(mmTOh[5])
DelThPlots.append(BIB[5])



fig, ax = plt.subplots()

h = ax.hist(ThPPlots,
        bins = 50,
        histtype = 'step',
        label = ['DiHiggs Processes',
                 '$\mu \mu \\rightarrow H^0 \mu \mu$',
                 '$\mu \mu \\rightarrow H^0 \\nu_m \\nu_m$',
                 '$\mu \mu \\rightarrow H^0$',
                 'BIB'])

plt.xlabel('$\\theta$ [rad]$')
plt.ylabel('Num Hits')
plt.yscale('log')
plt.title('$\\theta$ Distribution')

plt.savefig(f"ThDist.png")
