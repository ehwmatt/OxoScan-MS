# OxoScan-MS

A repository containing code used for analysis of glycoproteomics spectra.

## Functionality overview

The `glycoproteomics` Python library provides functions for working with and manipulating glycoproteomics spectra.
It was designed to be used within iPython notebooks, although can be used in other projects.

It is split into four parts:

- `io` for reading in spectra
- `spectrum` for binning, combining, and aligning spectra
- `peaks` for performing peak calling on the spectra, and extracting intensities from those peaks
- `plotting` for providing some basic plotting features within iPython notebooks

Spectra are read in and passed around between functions as hierarchical Python dictionaries, with the structure:

```
spectra_dict[rt_value][mz_value][ion_name] = intensity
```

### Example

This is some boilerplate code to import the various libraries and set up matplotlib:

```
import os
import glycoproteomics
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("png")
figure_size = (8, 4)
dpi = 80
```

Read in the spectrum and have a look to see which ions were quantified:

```
spectrum = glycoproteomics.io.read_spectrum_file("tests/data/spectrum.txt.gz")
ions = glycoproteomics.spectrum.list_ions(spectrum)
print(ions)
# ['138.055', '144.066', '168.066', '186.076', '204.087', '243.026', '274.092', '292.103', '308.098', '366.139', '405.079', '485.046', '512.197', '657.235']
```

In order to make the spectrum easier to work with, we bin the spectrum into RT and MZ bins: 

```
rt_x_bin_size = 0.025
mz_y_bin_size = 2.0

binned_spectrum = glycoproteomics.spectrum.bin(spectrum, rt_x_bin_size, mz_y_bin_size, np.mean)
```

We merge the spectra from individual ions into a matrix, which we can then plot:

```
ion_matrix, x_label, y_label = glycoproteomics.spectrum.to_matrix(binned_spectrum, ions)
glycoproteomics.plotting.plot_ion_matrix(ion_matrix, x_label, y_label, "spectrum.txt.gz", figure_size, dpi)
plt.show()
```

![Spectrum](/readme_images/spectrum.png)

Once in the matrix format, we can perform peak calling on the spectra to identify the top 10 peaks:

```
top_N_peaks = 10

# Peak quantification ellipse
x_radius = rt_x_bin_size * 3.0
y_radius = mz_y_bin_size * 5.0

# Peak exclusion ellipse (within which the centre of another peak will not be called)
x_radius_exclude = x_radius * 3.0
y_radius_exclude = y_radius * 2.0

peaks = glycoproteomics.peaks.find(ion_matrix, x_label, y_label, top_N_peaks, x_radius_exclude, x_radius_exclude)

glycoproteomics.plotting.plot_ion_matrix_with_peaks(
    ion_matrix, x_label, y_label, peaks, x_radius, y_radius, "spectrum.txt.gz - Top {} peaks".format(top_N_peaks), figure_size, dpi
)
plt.show()
```

![Spectrum with top 10 peaks](/readme_images/peaks.png)

## Getting set up

Set up a Python environment with `glycoproteomics` installed with the following commands.
This will then let you run the various workbooks which use the library.

```
virtualenv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install --upgrade -r requirements.txt -e .
```

## Running tests and profiling

Run the test suite with:

```
pytest -k "not profiling" --cov=glycoproteomics
```

To profile the code, and plot a graph of which functions take the most time, run:

```
pytest -k profiling
python -m gprof2dot -f pstats prof/peak_integration.out | dot -Tpdf -o prof/peak_integration.pdf
```

## Acknowledgements

Persistence code is cloned from a [repository](https://git.sthu.org/?p=persistence.git) by [Stefan Huber](https://www.sthu.org/code/codesnippets/imagepers.html). Licensed under version 3 of the GNU Lesser General Public License.
