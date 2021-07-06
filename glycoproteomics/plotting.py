import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def plot_ion_matrix(ion_matrix, x_label, y_label, title, figsize, dpi):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    c = ax.pcolormesh(
        np.rot90(np.array([x_label] * len(y_label)), 3),
        np.array([y_label] * len(x_label)),
        np.log10(ion_matrix + 1.0),
        cmap="viridis",
        shading="nearest",
    )
    ax.set_title(title)
    ax.set_xlabel("RT")
    ax.set_ylabel("MZ")
    fig.colorbar(c, ax=ax)
    return fig, ax


def plot_ion_matrix_with_peaks(
    ion_matrix, x_label, y_label, peaks, x_radius, y_radius, title, figsize, dpi
):
    fig, ax = plot_ion_matrix(ion_matrix, x_label, y_label, title, figsize, dpi)
    for idx in range(len(peaks)):
        ax.add_patch(
            Ellipse(
                peaks[idx][0],
                x_radius * 2.0,
                y_radius * 2.0,
                facecolor="none",
                edgecolor="red",
            )
        )
        ax.text(
            peaks[idx][0][0] - x_radius,
            peaks[idx][0][1],
            str(idx + 1),
            horizontalalignment="right",
            verticalalignment="center",
            fontsize=6,
            color="white",
        )
    return fig, ax
