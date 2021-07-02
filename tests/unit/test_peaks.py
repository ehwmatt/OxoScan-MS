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
