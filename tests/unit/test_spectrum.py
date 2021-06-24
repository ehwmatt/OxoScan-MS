import glycoproteomics
import numpy as np
import os


def test_bin_spectrum():
    spectrum = {
        0.00350036: {
            780.8955: {
                "204.087": 1.0,
                "186.076": 2.0,
                "168.066": 3.0,
            }
        },
        0.00450542: {
            782.887: {
                "204.087": 4.0,
                "186.076": 5.0,
                "168.066": 6.0,
            }
        },
        0.0049: {
            783.887: {
                "204.087": 7.0,
                "186.076": 8.0,
                "168.066": 9.0,
            }
        },
    }
    assert glycoproteomics.spectrum.bin(spectrum, 0.02, 4.0, np.mean) == {
        0.01: {
            782.0: {
                "204.087": 4.0,
                "186.076": 5.0,
                "168.066": 6.0,
            }
        }
    }
    assert glycoproteomics.spectrum.bin(spectrum, 0.02, 4.0, np.max) == {
        0.01: {
            782.0: {
                "204.087": 7.0,
                "186.076": 8.0,
                "168.066": 9.0,
            }
        }
    }
    assert glycoproteomics.spectrum.bin(spectrum, 0.002, 2.0, np.mean) == {
        0.003: {
            781.0: {
                "204.087": 1.0,
                "186.076": 2.0,
                "168.066": 3.0,
            },
        },
        0.005: {
            783.0: {
                "204.087": 5.5,
                "186.076": 6.5,
                "168.066": 7.5,
            }
        },
    }


def test_bin_test_spectrum(rootdir):
    spectrum_path = os.path.join(rootdir, "tests", "data", "spectrum.txt.gz")
    spectrum = glycoproteomics.io.read_spectrum_file(spectrum_path)
    binned_spectrum = glycoproteomics.spectrum.bin(spectrum, 1.0, 100.0, np.mean)
    assert list(binned_spectrum.keys()) == [
        0.5,
        1.5,
        2.5,
        3.5,
        4.5,
        5.5,
        6.5,
        7.5,
        8.5,
        9.5,
        10.5,
        11.5,
        12.5,
        13.5,
        14.5,
    ]
    assert list(binned_spectrum[0.5].keys()) == [
        750.0,
        850.0,
        950.0,
        1050.0,
        1150.0,
        1250.0,
        1350.0,
        1450.0,
        1550.0,
        1650.0,
    ]
    assert binned_spectrum[0.5][750.0] == {
        "204.087": 0.0,
        "186.076": 0.0,
        "168.066": 0.0,
        "366.139": 0.0,
        "144.066": 0.0,
        "138.055": 0.0,
        "512.197": 0.007032763636363635,
        "292.103": 0.0,
        "274.092": 0.0,
        "657.235": 0.0,
        "243.026": 0.0,
        "405.079": 0.0,
        "485.046": 0.0,
        "308.098": 0.0,
    }


def test_combine_spectra():
    multiple_spectra = {
        "spectrum1.txt": {
            0.05: {
                780.0: {
                    "204.087": 1.0,
                    "186.076": 2.0,
                    "168.066": 3.0,
                }
            },
            0.07: {
                782.0: {
                    "204.087": 4.0,
                    "186.076": 5.0,
                    "168.066": 6.0,
                }
            },
        },
        "spectrum2.txt": {
            0.05: {
                780.0: {
                    "204.087": 1.0,
                    "186.076": 2.0,
                    "168.066": 3.0,
                }
            },
            0.06: {
                782.0: {
                    "204.087": 4.0,
                    "186.076": 5.0,
                    "168.066": 6.0,
                }
            },
        },
        "spectrum3.txt": {
            0.05: {
                780.0: {
                    "204.087": 1.0,
                    "186.076": 2.0,
                    "168.066": 3.0,
                }
            },
            0.07: {
                782.0: {
                    "204.087": 4.0,
                    "186.076": 5.0,
                    "168.066": 6.0,
                }
            },
        },
    }
    assert glycoproteomics.spectrum.combine(multiple_spectra, np.sum) == {
        0.05: {
            780.0: {
                "204.087": 3.0,
                "186.076": 6.0,
                "168.066": 9.0,
            }
        },
        0.06: {
            782.0: {
                "204.087": 4.0,
                "186.076": 5.0,
                "168.066": 6.0,
            }
        },
        0.07: {
            782.0: {
                "204.087": 8.0,
                "186.076": 10.0,
                "168.066": 12.0,
            }
        },
    }


def test_spectrum_to_matrix():
    spectrum = {
        0.00350036: {
            780.8955: {
                "204.087": 1.0,
                "186.076": 2.0,
                "168.066": 3.0,
            }
        },
        0.00450542: {
            782.887: {
                "204.087": 4.0,
                "186.076": 5.0,
                "168.066": 6.0,
            }
        },
        0.0049: {
            783.887: {
                "204.087": 7.0,
                "186.076": 8.0,
                "168.066": 9.0,
            }
        },
    }
    result = glycoproteomics.spectrum.to_matrix(spectrum, "204.087")
    np.testing.assert_array_equal(
        result[0], np.array([[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 7.0]])
    )
    assert result[1] == [0.00350036, 0.00450542, 0.0049]
    assert result[2] == [780.8955, 782.887, 783.887]
