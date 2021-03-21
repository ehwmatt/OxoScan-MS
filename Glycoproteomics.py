python_env_path = '/camp/lab/ralserm/working/Glycoproteomics/Python Environment/' # this is to load tqdm which isnt installed on CAMP

# Load parameters setup in ipynb by running setup_parameters.py or by default in setup_parameters.py

from setup_parameters import *

print('Imported updated parameters from setup_parameters.py into Glycoproteomics.py')
# for loading massspec reports

import collections
import csv

# for preprocessing

import math
import numpy as np
from matplotlib import cm

# for persistance

import warnings
from multiprocessing import Pool
from os import listdir
from os.path import join, isfile
import numpy as np
from python_env_path import tqdm

# For visualization

import numpy as np
import matplotlib.pyplot as plt

#####################################################################
### Functions for loading extracted mass spec reports from DIA-NN ###
#####################################################################

def load(filename):
    dict_data = {}
    ion_labels = []
    with open(filename, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        row = next(reader)
        for i in range(3, len(row)):
            ion_labels.append(row[i])
        for row in reader:
            rt = float(row[0])
            mz = (float(row[1]) + float(row[2])) / 2
            ion_dict = {}
            for i in range(3, len(row)):
                ion_dict[ion_labels[i - 3]] = float(row[i])
            if rt not in dict_data:
                dict_data[rt] = {}
            row = dict_data[rt]
            row[mz] = ion_dict

    return dict_data


def get_all_ions(dict_data):
    all_ions = []
    for rt in dict_data:
        rt_dict = dict_data[rt]
        for mz in rt_dict:
            ions = rt_dict[mz]
            if isinstance(ions, collections.Mapping):
                for ion in ions:
                    all_ions.append(ion)
            return all_ions


def save(filename, dict_data):
    column_headers = ['RT', 'MZ']
    all_ions = get_all_ions(dict_data)
    for ion in all_ions:
        column_headers.append(ion)
    with open(filename, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(column_headers)
        for rt in dict_data:
            rt_dict = dict_data[rt]
            for mz in rt_dict:
                ions = rt_dict[mz]
                row = [rt, mz]
                if isinstance(ions, collections.Mapping):
                    for ion in all_ions:
                        row.append(ions[ion])
                else:
                    row.append(ions)
                writer.writerow(row)

############################################################################
### Functions for preprocessing mass spec reports loaded as dictionaries ###
############################################################################

# Function to get dimensions of rt vs mz from a mass spectra dictionary

def get_data_size(dict_data):
    width = len(dict_data)
    height = 0
    for rt in dict_data:
        rt_dict = dict_data[rt]
        height = len(rt_dict)
    return height, width


# Function to convert loaded ms spectra dictionary to 3 lists x (retention time) , y(mz), intensity

def convert_to_xy(dict_data):
    x_data = []
    y_data = []
    intensity_data = []
    for rt in dict_data:
        rt_dict = dict_data[rt]
        for mz in rt_dict:
            ions = rt_dict[mz]
            intensity = 0
            for ion in ions:
                intensity += ions[ion]
            intensity /= len(ions)
            x_data.append(rt)
            y_data.append(mz)
            intensity_data.append(intensity)
    return x_data, y_data, intensity_data

# Function to convert loaded ms spectra dictionary to 3 lists x (retention time) , y(mz), colour scale

def convert_to_xyc(dict_data, select_ions=None):
    x_data = []
    y_data = []
    color_data = []
    for rt in dict_data:
        rt_dict = dict_data[rt]
        for mz in rt_dict:
            ions = rt_dict[mz]
            if select_ions:
                loop_ions = select_ions
            else:
                loop_ions = ions
            color_labels = cm.rainbow(np.linspace(0, 1, len(loop_ions)))
            color = np.zeros_like(color_labels[0])
            for i, ion in enumerate(loop_ions):
                intensity = ions[ion]
                if intensity:
                    intensity = 1 / (2 - math.log10(min(intensity, 10)))
                color += (1 - (1 - color_labels[int(i)]) * intensity)
            color /= len(loop_ions)
            x_data.append(rt)
            y_data.append(mz)
            color_data.append(color)
    return x_data, y_data, color_data


def convert_to_matrix(dict_data, ion):
    x_data = []
    y_data = []
    for rt in dict_data:
        x_data.append(rt)
        rt_dict = dict_data[rt]
        for mz in rt_dict:
            y_data.append(mz)
    sorted_x_data = sorted(list(set(x_data)))
    sorted_y_data = sorted(list(set(y_data)))
    ion_data = np.zeros((len(sorted_x_data), len(sorted_y_data)))
    for rt in dict_data:
        rt_dict = dict_data[rt]
        for mz in rt_dict:
            ion_data[sorted_x_data.index(rt), sorted_y_data.index(mz)] = rt_dict[mz][ion]
    return sorted_x_data, sorted_y_data, ion_data

# Function to resample ( bin spectra into matrix ) according to parameters specified in setup_parameters.py )

def resample_to_matrix(dict_data, min_rt, max_rt, step_rt, min_mz, max_mz, step_mz, ion, normalise=False):
    w = int((max_rt - min_rt) / step_rt)
    h = int((max_mz - min_mz) / step_mz)
    matrix = np.zeros((w, h))
    for rt in dict_data:
        rt_dict = dict_data[rt]
        for mz in rt_dict:
            x = int((rt - min_rt) / step_rt)
            y = int((mz - min_mz) / step_mz)
            if x < 0:
                x = 0
            if x >= w:
                x = w - 1
            if y < 0:
                y = 0
            if y >= h:
                y = h - 1
            matrix[x, y] += rt_dict[mz][selected_ion]

    if normalise:
        matrix /= np.mean(matrix)

    return matrix

###########################################
### Functions for persistance algorithm ###
###########################################
"""UnionFind.py

Union-find data structure. Based on Josiah Carlson's code,
http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
with significant additional changes by D. Eppstein.
"""


class UnionFind:

    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def add(self, object, weight):
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = weight

    def __contains__(self, object):
        return object in self.parents

    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            assert(False)
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure.

        """
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r], r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.parents[r] = heaviest

# Source: https://www.sthu.org/code/codesnippets/imagepers.html

"""A simple implementation of persistent homology on 2D images."""

__author__ = "Stefan Huber <shuber@sthu.org>"

def get(im, p):
    return im[p[0]][p[1]]


def iter_neighbors(p, w, h):
    y, x = p

    # 8-neighborship
    neigh = [(y+j, x+i) for i in [-1, 0, 1] for j in [-1, 0, 1]]
    # 4-neighborship
    # neigh = [(y-1, x), (y+1, x), (y, x-1), (y, x+1)]

    for j, i in neigh:
        if j < 0 or j >= h:
            continue
        if i < 0 or i >= w:
            continue
        if j == y and i == x:
            continue
        yield j, i


def persistence(im):
    h, w = im.shape

    # Get indices orderd by value from high to low
    indices = [(i, j) for i in range(h) for j in range(w)]
    indices.sort(key=lambda p: get(im, p), reverse=True)

    # Maintains the growing sets
    uf = UnionFind()

    groups0 = {}

    def get_comp_birth(p):
        return get(im, uf[p])

    # Process pixels from high to low
    for i, p in enumerate(indices):
        v = get(im, p)
        ni = [uf[q] for q in iter_neighbors(p, w, h) if q in uf]
        nc = sorted([(get_comp_birth(q), q) for q in set(ni)], reverse=True)

        if i == 0:
            groups0[p] = (v, v, None)

        uf.add(p, -i)

        if len(nc) > 0:
            oldp = nc[0][1]
            uf.union(oldp, p)

            # Merge all others with oldp
            for bl, q in nc[1:]:
                if uf[q] not in groups0:
                    #print(i, ": Merge", uf[q], "with", oldp, "via", p)
                    groups0[uf[q]] = (bl, bl-v, p)
                uf.union(oldp, q)

    groups0 = [(k, groups0[k][0], groups0[k][1], groups0[k][2]) for k in groups0]
    groups0.sort(key=lambda g: g[2], reverse=True)

    return groups0

################################################################################################
### Functions for finding rt-mz peaks (with persistance algorithm) & estimating peak volumes ###
################################################################################################

def get_sigma(matrix, pos, delta):
    x = pos[0] + delta[0]
    y = pos[1] + delta[1]
    w = matrix.shape[0]
    h = matrix.shape[1]
    if 0 <= x < w and 0 <= y < h:
        value0 = matrix[pos]
        value = matrix[x, y]
        if 0 < value < value0:
            sigma = np.linalg.norm(delta) / np.sqrt(np.log(value / value0) / -0.5)
            return sigma
    return np.nan


def get_integrated(matrix, pos, delta):
    area = 0
    w = matrix.shape[0]
    h = matrix.shape[1]
    mask = [(pos[0] + i, pos[1] + j) for i in range(-delta, delta + 1) for j in range(-delta, delta + 1)]
    for x, y in mask:
        if 0 <= x < w and 0 <= y < h:
            area += matrix[x, y]
    return area


def estimate_integrated(matrix, pos):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sigmax = np.nanmean([get_sigma(matrix, pos, (-1, 0)), get_sigma(matrix, pos, (1, 0))])
        sigmay = np.nanmean([get_sigma(matrix, pos, (0, -1)), get_sigma(matrix, pos, (0, 1))])

    if np.isnan(sigmax):
        sigmax = 0
    if np.isnan(sigmay):
        sigmay = 0

    if not sigmax and not sigmay:
        area = matrix[pos]
    elif not sigmax:
        area = sigmay * np.sqrt(2 * np.pi) * matrix[pos]
    elif not sigmay:
        area = sigmax * np.sqrt(2 * np.pi) * matrix[pos]
    else:
        area = sigmax * sigmay * 2 * np.pi * matrix[pos]

    return area


def calc_peaks(filename):
    data = load(filename)
    matrix = resample_to_matrix(data,
                                min_rt=min_rt, max_rt=max_rt, step_rt=step_rt,
                                min_mz=min_mz, max_mz=max_mz, step_mz=step_mz,
                                ion=selected_ion, normalise=True)
    peaks0 = persistence(matrix)
    peaks = []
    # take [npeaks] top number of peaks
    for peak0 in peaks0[0:npeaks]:
        x = (peak0[0][0] + 0.5) * step_rt + min_rt
        y = (peak0[0][1] + 0.5) * step_mz + min_mz
        intensity = estimate_integrated(matrix, peak0[0])
        peaks.append((x, y, intensity))
    return np.asarray(peaks)


def process_peaks(filename):
    peaks = calc_peaks(filename)
    outfile = filename.replace(file_extension_ms, file_extension_peaks)
    np.savetxt(outfile, peaks, fmt='%f', delimiter=',')
    return peaks

######################################
### Functions for Visualizing data ###
######################################

def plot_init():
    plt.figure(dpi=600)
    plt.xlim(min_rt, max_rt)
    plt.ylim(min_mz, max_mz)


def plot_peaks(peaks, peak_index=2, base_size=20):
    x_data = []
    y_data = []
    intensity = []
    for peak in peaks:
        x_data.append(peak[0])
        y_data.append(peak[1])
        intensity.append(peak[peak_index])
    plt.scatter(x_data, y_data, s=base_size * np.sqrt(intensity), alpha=0.5, marker='.', edgecolor='')

###################################
### Functions for merging peaks ###
###################################


# ... coming soon ...

