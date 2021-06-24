import glycoproteomics.io
import os


def test_read_spectrum_file(tmp_path):
    test_spectrum = os.path.join(tmp_path, "spectrum.txt")
    with open(test_spectrum, "w") as spectrum:
        spectrum.write("RT\tWindow.Low\tWindow.High\t204.087\t186.076\t168.066\n")
        spectrum.write("0.00450036\t779.9\t781.891\t1\t2\t3\n")
        spectrum.write("0.00450542\t781.891\t783.883\t4\t5\t6\n")
    assert glycoproteomics.io.read_spectrum_file(test_spectrum) == {
        0.00450036: {
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
    }


def test_read_gzipped_spectrum_file(rootdir):
    test_spectrum = os.path.join(rootdir, "tests", "data", "spectrum.txt.gz")
    assert isinstance(glycoproteomics.io.read_spectrum_file(test_spectrum), dict)


def test_read_spectra_directory(tmp_path):
    for idx in range(3):
        test_spectrum = os.path.join(tmp_path, "spectrum{}.txt".format(idx + 1))
        with open(test_spectrum, "w") as spectrum:
            spectrum.write("RT\tWindow.Low\tWindow.High\t204.087\t186.076\t168.066\n")
            spectrum.write("0.00450036\t779.9\t781.891\t1\t2\t3\n")
            spectrum.write("0.00450542\t781.891\t783.883\t4\t5\t6\n")
    assert glycoproteomics.io.read_spectra_directory(tmp_path) == {
        "spectrum1.txt": {
            0.00450036: {
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
        },
        "spectrum2.txt": {
            0.00450036: {
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
        },
        "spectrum3.txt": {
            0.00450036: {
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
        },
    }
