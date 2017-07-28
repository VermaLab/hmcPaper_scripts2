import peakProfiles as pp
import numpy as np
import os
import sys

reload(pp)

resFolder = "hpneMcSites_hpneAtac"
resFolder = 'hpne_mc_100_hpneAtac'
resFolder = 'jd_mc_100_jdAtac'
resFolder = sys.argv[1]

#peakProfs = np.load(os.path.join(resFolder, "peakProfiles.npy"))
prof = np.load(os.path.join(resFolder, "summedPP.npy"))
smallPeaks = np.load(os.path.join(resFolder, "binnedPP.npy"))

pp.savePeakPlots(os.path.join(resFolder, "colocPlot.png"), prof, smallPeaks)
pp.savePeakPlots_part(os.path.join(resFolder, "colocPlot_20.png"), prof, smallPeaks, throwPortion = .20)
pp.savePeakPlots_part(os.path.join(resFolder, "colocPlot_50.png"), prof, smallPeaks, throwPortion = .50)


