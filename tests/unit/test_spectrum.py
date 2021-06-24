import glycoproteomics
import numpy as np


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
