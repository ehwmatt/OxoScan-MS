import glycoproteomics
import numpy as np


def test_peaks_find():
    ion_matrix = np.array(
        [
            [1, 1, 1, 1, 1],
            [1, 2, 2, 2, 1],
            [1, 2, 3, 2, 1],
            [1, 2, 2, 2, 4],
        ]
    )
    x_labels = [0.02, 0.04, 0.06, 0.08]
    y_labels = [800.0, 802.0, 804.0, 806.0, 808.0]
    # Test basic calling
    assert glycoproteomics.peaks.find(
        ion_matrix, x_labels, y_labels, 10, 0.02, 2.0
    ) == [((0.08, 808.0), 4, 4), ((0.06, 804.0), 3, 1)]
    # Test calling where peaks which are too close are excluded
    assert glycoproteomics.peaks.find(
        ion_matrix, x_labels, y_labels, 10, 0.06, 6.0
    ) == [((0.08, 808.0), 4, 4)]
    assert glycoproteomics.peaks.find(
        ion_matrix, x_labels, y_labels, 10, 0.01, 4.0
    ) == [((0.08, 808.0), 4, 4), ((0.06, 804.0), 3, 1)]
    # Test calling fewer peaks
    assert glycoproteomics.peaks.find(ion_matrix, x_labels, y_labels, 1, 0.02, 2.0) == [
        ((0.08, 808.0), 4, 4)
    ]


def test_peaks_integrate():
    ion_matrix = np.array(
        [
            [1, 1, 1, 1, 1],
            [1, 2, 2, 2, 1],
            [1, 2, 3, 2, 1],
            [1, 2, 2, 2, 4],
        ]
    )
    x_labels = [0.2, 0.4, 0.6, 0.8]
    y_labels = [800.0, 802.0, 804.0, 806.0, 808.0]
    peaks = [((0.8, 808.0), 4, 4), ((0.6, 804.0), 3, 1)]
    assert glycoproteomics.peaks.integrate(
        ion_matrix, x_labels, y_labels, peaks, 0.3, 3.0, np.max
    ) == [4.0, 3.0]
    # Test different integration function
    assert glycoproteomics.peaks.integrate(
        ion_matrix, x_labels, y_labels, peaks, 0.3, 3.0, np.sum
    ) == [9.0, 19.0]
    # Test different radii
    assert glycoproteomics.peaks.integrate(
        ion_matrix, x_labels, y_labels, peaks, 0.3, 5.0, np.sum
    ) == [11.0, 21.0]


def test_peaks_rt_move():
    peaks = [((0.8, 808.0), 4, 4), ((0.6, 804.0), 3, 1)]
    rt_alignment = {0.8: 0.7, 0.6: 0.5}
    assert glycoproteomics.peaks.rt_move(peaks, rt_alignment) == [
        ((0.7, 808.0), 4, 4),
        ((0.5, 804.0), 3, 1),
    ]
