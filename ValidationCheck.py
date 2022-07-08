import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle

fig, (ax1, ax2) = plt.subplots(ncols=2)


with open("Expected_Hists.pickle", 'rb') as f:
    Exp_Hists = pickle.load(f)

with open("New_Hists.pickle", 'rb') as f:
        New_Hists = pickle.load(f)

TimePlots = []
TimePlots.append(Exp_Hists[0])
TimePlots.append(New_Hists[0])

ThetaPlots = []
ThetaPlots.append(Exp_Hists[1])
ThetaPlots.append(New_Hists[1])

ns, bins, patches = ax1.hist(TimePlots,
                            histtype = 'stepfilled',
                            bins = 30,
                            alpha = 0.2,
                            label=['Expected','New'])

ax1.legend()

ns, bins, patches = ax2.hist(ThetaPlots,
                                histtype = 'stepfilled',
                                bins = 30,
                                alpha = 0.2,
                                label=['Expected','New'])

plt.savefig(f"Validation_Hists")
