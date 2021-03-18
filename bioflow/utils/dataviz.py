"""
All credits for the implementation and suggestions go to sega_sai (stackoverflow):
http://stackoverflow.com/questions/10439961/efficiently-create-a-density-plot-for-high-density-regions-points-for-sparse-re

As a datavisualization routine, this is a little bit hard to unitttest. As such, it has been
excluded from the unittesting framework
"""
# from scipy.sparse import lil_matrix, triu
# from bioflow.utils.linalg_routines import normalize_laplacian
import matplotlib.pyplot as plt
import numpy as np
from scipy import histogram2d
from scipy.stats import gaussian_kde
import os
from bioflow.configs.main_configs import output_location
from typing import Any, Union, TypeVar, NewType, Tuple, List


def better_2d_density_plot(x_data, y_data, threshold=3, bins=(100, 100)):
    """
    Takes x and y coordinates of the data points and creates a 2D density plot. Everything below
    threshold is represented as individual points, bins tuple is the number of bins along each axis.
    Essentially is a wrapper around scipy's histogram2d. Does not show, so a pyplot.show() is
    required to show the result of plotting

    :param x_data: x coordinates of the data points
    :param y_data: y coordinates of the data points
    :param threshold: if there are less than that number of points in a bin, points are
    represented individually
    :param bins: number of bins along each axis
    """
    xy_range = [[min(x_data), max(x_data)], [min(y_data), max(y_data)]]
    distortion = (xy_range[1][1] - xy_range[1][0]) / \
        (xy_range[0][1] - xy_range[0][0])
    x_data = x_data * distortion

    xy_range = [[min(x_data), max(x_data)], [min(y_data), max(y_data)]]
    hh, loc_x, loc_y = histogram2d(x_data, y_data, range=xy_range, bins=bins)
    pos_x = np.digitize(x_data, loc_x)
    pos_y = np.digitize(y_data, loc_y)

    ind = (pos_x > 0) & (pos_x <= bins[0]) & (pos_y > 0) & (pos_y <= bins[1])
    # values of the histogram where the points are
    hh_sub = hh[pos_x[ind] - 1, pos_y[ind] - 1]
    x_dat1 = x_data[ind][hh_sub < threshold]  # low density points
    y_dat1 = y_data[ind][hh_sub < threshold]
    hh[hh < threshold] = np.nan  # fill the areas with low density by NaNs

    plt.imshow(np.flipud(hh.T), cmap='jet',
               extent=np.array(xy_range).flatten(), interpolation='none')
    plt.plot(x_dat1, y_dat1, '.')


def violin_plot(axis, data_, position_, box_plot=False):
    """
    Creates a violin plot along the axis in the figure

    :param axis: axis from an existing figure where violin plot should be created
    :param data_: data used to create a violin plot
    :param position_: tuple of positions limiting the violin plot and boxplot
    :param box_plot: if true, overlays a boxplot atop the violin plot
    :return:
    """
    dist = max(position_) - min(position_)
    w = min(0.15 * max(dist, 1.0), 0.5)
    for d, p in zip(data_, position_):
        kernel_density = gaussian_kde(d)
        low_bound = kernel_density.dataset.min()
        upper_bound = kernel_density.dataset.max()
        violing_support = np.arange(
            low_bound,
            upper_bound,
            (upper_bound - low_bound) / 100.)
        violin_profile = kernel_density.evaluate(violing_support)
        violin_profile = violin_profile / violin_profile.max() * w
        # scaling the violin to the available space
        axis.fill_betweenx(
            violing_support,
            p,
            violin_profile + p,
            facecolor='y',
            alpha=0.3)
        axis.fill_betweenx(
            violing_support, p, -violin_profile + p, facecolor='y', alpha=0.3)
    if box_plot:
        axis.boxplot(data_, notch=1, positions=position_, vert=1)


def kde_compute(bi_array, bin_no=30, samples=10, show=True):
    """
    Computes a kernel density estimator and plots it in case show is true.

    :param bi_array: array of points coordinates that are used to compute the kde
    :param bin_no: number of bins used to
    :param samples: number of samples used to estimate the kde
    :param show: if kernel density estimator should be plotted or not. If yes, an additional
    plt.show() is required to show the result after that method
    :return:
    """
    repeated_sample_correction = bi_array.shape[1] / float(samples)
    x, y = bi_array

    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data
    # extents
    k = gaussian_kde(bi_array)
    xi, yi = np.mgrid[x.min():x.max():bin_no * 1j, y.min():y.max():bin_no * 1j]
    zi = np.tanh(k(np.vstack([xi.flatten(), yi.flatten()]))
                 * repeated_sample_correction)

    if show:
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')

    return lambda x_: np.tanh(k(x_) * repeated_sample_correction)


# def view_laplacian_off_terms(non_normalized_laplacian):
#     """
#     Shows a log-10 histogram of distribution of off-diagonal terms
#
#     :param non_normalized_laplacian:
#     :return:
#     """
#     # if we revive, the line below will have to module that will be revived
#     normalized_laplacian = normalize_laplacian(non_normalized_laplacian)
#     triangular_upper = lil_matrix(triu(normalized_laplacian))
#     triangular_upper.setdiag(0)
#     pre_arr = -triangular_upper[triangular_upper.nonzero()].toarray().flatten()
#     arr = np.log10(pre_arr)
#     plt.hist(arr, bins=100, log=True, histtype='step')
#     plt.show()

# TRACING: [run path] pipe hdd save destination here (0)
def render_2d_matrix(matrix: np.array,
                     name: str,
                     destination: str = '') -> None:
    """
    Subroutine required by the rendering wrapper.

    :param matrix:
    :param name:
    :return:
    """
    plt.title(name)
    plt.imshow(matrix, interpolation='nearest')
    plt.colorbar()
    # TRACING: [run path] here is where we save the clustering
    plt.savefig(os.path.join(output_location, name + '.png'))


if __name__ == "__main__":
    from numpy.random import normal
    N = 1e5
    x_dat, y_dat = np.random.normal(size=N), np.random.normal(1, 0.6, size=N)
    better_2d_density_plot(x_dat, y_dat)
    plt.show()

    pos = list(range(5))
    data = [normal(size=100) for i in pos]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    violin_plot(ax, data, pos, box_plot=1)
    plt.show()

    np.random.seed(1977)

    # Generate 200 correlated x,y points
    data = np.random.multivariate_normal([0, 0], [[1, 0.5], [0.5, 3]], 200)
    kde_compute(data.T, bin_no=20)
    plt.show()
