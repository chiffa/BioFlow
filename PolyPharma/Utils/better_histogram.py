"""
All credits for the implementation and suggestions go to sega_sai (stackoverflow):
http://stackoverflow.com/questions/10439961/efficiently-create-a-density-plot-for-high-density-regions-points-for-sparse-re
"""
__author__ = 'ank'


import matplotlib.pyplot as plt
import numpy as np
from scipy import histogram2d


def better2D_desisty_plot(xdat, ydat, thresh=3, bins = (100, 100)):
    xyrange = [[min(xdat), max(xdat)], [min(ydat), max(ydat)]]
    distortion = (xyrange[1][1]-xyrange[1][0])/(xyrange[0][1]-xyrange[0][0])
    xdat = xdat * distortion

    xyrange = [[min(xdat), max(xdat)], [min(ydat), max(ydat)]]
    hh, locx, locy = histogram2d(xdat, ydat, range=xyrange, bins=bins)
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)

    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh] # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan # fill the areas with low density by NaNs

    plt.imshow(np.flipud(hh.T), cmap='jet', extent=np.array(xyrange).flatten(), interpolation='none')
    plt.plot(xdat1, ydat1, '.')


if __name__ == "__main__":
    N = 1e5
    xdat, ydat = np.random.normal(size=N), np.random.normal(1, 0.6, size=N)
    better2D_desisty_plot(xdat, ydat)
    plt.show()