# Glycoproteomics

A repository containing code used for analysis of glycoproteomics spectra.

## Getting set up

Set up a Python environment with `glycoproteomics` installed with the following commands.
This will then let you run the various workbooks which use the library.

```
virtualenv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install --upgrade -r requirements.txt -e .
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
