import numpy as np
import matplotlib as mpl
from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
from pylab import *
import seaborn as sns
sns.set(style="whitegrid")
sns.set_context("notebook", font_scale=2.0, rc={"lines.linewidth":1.0})
sns.axes_style(rc={"axes.edgecolor":'black'})
matplotlib.rc('font', size=17)
matplotlib.rc('figure', autolayout=True)
matplotlib.rc('axes', edgecolor='0.57')
matplotlib.rc('grid', color='0.57')
import healpy as hp
import seaborn as sns

cmb = hp.read_map("/home/mabitbol/mapmaking/simulator/planck_data/cmb_map512.fits", field=(0,1,2))
norm1 = 1e6
Tcmb = cmb[0]*norm1
Qcmb = cmb[1]*norm1
Ucmb = cmb[2]*norm1
dir1 = "/home/mabitbol/mapmaking/simulator/simdata/white25/"
dir2 = "/home/mabitbol/mapmaking/simulator/simdata/red25/"
f1 = np.load(dir1+"in_maps.npz")
f2 = np.load(dir2+"in_maps.npz")
f1.keys()
I1 = f1['I']
I2 = f2['I']
hits1 = f1['hits']
mask = hits1 < 1
Tcmb[mask] = hp.UNSEEN
des1 = hp.read_map(dir1+"map_o400.fits")
destriped = hp.read_map(dir2+"map_o40.fits")
naive1 = hp.read_map(dir1+"naive_map_o400.fits")
naive = hp.read_map(dir2+"naive_map_o40.fits")

figure()
hp.mollview(Tcmb, min=-400, max=400)
title('CMB')
savefig('cmb')

figure()
hp.mollview(naive, min=-1000, max=1000)
title('Binned Map')
savefig('stripedmap')

mask = destriped == hp.UNSEEN
rmean = np.mean(destriped[~mask])
destriped -= rmean
figure()
hp.mollview(destriped, min=-400, max=400)
title('Destriped')
savefig('destriped')

bestcase = naive1 - Tcmb
mask = naive1 == hp.UNSEEN
mask |= Tcmb == hp.UNSEEN
bestcase[mask] = hp.UNSEEN

cut = 1000
residuals = destriped - Tcmb
mask = Tcmb == hp.UNSEEN
mask |= destriped == hp.UNSEEN
residuals[mask] = hp.UNSEEN
print "mean offset: ", rmean
mask = np.abs(residuals) > cut
print "residual standard deviation: ", np.std(residuals[~mask])
figure()
sns.distplot(residuals[~mask], kde=False, label="destriped-cmb", color='g')
mask = np.abs(bestcase) > cut
sns.distplot(bestcase[~mask], kde=False, label="white-cmb", color='r')
print "best case standard deviation: ", np.std(bestcase[~mask])
mask = naive == hp.UNSEEN
mask |= Tcmb == hp.UNSEEN
resnaive = naive - Tcmb
resnaive[mask] = hp.UNSEEN
mask = np.abs(resnaive) > 1000
sns.distplot(resnaive[~mask], kde=False, label='binned-cmb', color='b')
print "worst case standard deviation: ", np.std(resnaive[~mask])
xlabel('$\mu K$')
ylabel('number')
title('Histogram Residuals to CMB')
legend()
savefig('residualhist')
