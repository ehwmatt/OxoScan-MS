# Glycoproteomics

A repository containing code used for analysis of glycoproteomics spectra.

## Installing

Install Python requirements with:

```
pip install -r requirements.txt
```

## Running tests

Run the test suite with:

```
pytest -k "not profiling" --cov=glycoproteomics
```

## Development

To profile the code, and plot a graph of which functions take the most time, run:

```
pytest -k profiling
python -m gprof2dot -f pstats prof/peak_integration.out | dot -Tpdf -o prof/peak_integration.pdf
```

## Acknowledgements

Persistence code is cloned from a [repository](https://git.sthu.org/?p=persistence.git) by [Stefan Huber](https://www.sthu.org/code/codesnippets/imagepers.html). Licensed under version 3 of the GNU Lesser General Public License.
