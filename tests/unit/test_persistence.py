import numpy as np
import glycoproteomics.persistence.imagepers


def test_single_peak():
    assert (
        glycoproteomics.persistence.imagepers.persistence(
            np.array(
                [
                    [1, 1, 1, 1, 1],
                    [1, 2, 2, 2, 1],
                    [1, 2, 3, 2, 1],
                    [1, 2, 2, 2, 1],
                    [1, 1, 1, 1, 1],
                ]
            )
        )
        == [((2, 2), 3, 3, None)]
    )


def test_multi_peak():
    assert (
        glycoproteomics.persistence.imagepers.persistence(
            np.array(
                [
                    [1, 1, 1, 1, 1],
                    [1, 2, 2, 2, 1],
                    [1, 2, 3, 2, 1],
                    [1, 2, 2, 2, 1],
                    [1, 1, 1, 1, 4],
                ]
            )
        )
        == [((4, 4), 4, 4, None), ((2, 2), 3, 1, (3, 3))]
    )
