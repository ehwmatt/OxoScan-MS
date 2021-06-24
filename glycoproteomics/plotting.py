import matplotlib.pyplot as plt
import numpy as np


def plot_ion_matrix(ion_matrix):
    plt.imshow(np.rot90(np.log10(ion_matrix + 1.0)))
