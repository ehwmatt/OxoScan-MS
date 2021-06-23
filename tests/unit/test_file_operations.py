import glycoproteomics
import os


def test_load(tmp_path):
    test_spectra = os.path.join(tmp_path, "spectra.txt")
    with open(test_spectra, "w") as spectra:
        spectra.write("RT\tWindow.Low\tWindow.High\t204.087\t186.076\t168.066\n")
        spectra.write("0.00450036\t779.9\t781.891\t1\t2\t3\n")
        spectra.write("0.00450542\t781.891\t783.883\t4\t5\t6\n")
    assert glycoproteomics.read_spectra_file(test_spectra) == {
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
