# Simulating and analyzing a theoretical wine grape demography

This project is setup to perform large scale population genetic simulations and perform time-series analysis of the data.

## Clonal package
clonal is a [pybind11](https://pybind11.readthedocs.io/en/stable/) extension module for [fwdpy11](https://github.com/molpopgen/fwdpy11)
### clonal dependencies

* [Python3](https://www.python.org/download/releases/3.0/)
* [pybind11](https://pybind11.readthedocs.io/en/stable/) 
* [fwdpy11](https://github.com/molpopgen/fwdpy11)
* [pylibseq](https://github.com/molpopgen/pylibseq)
* [numpy](http://www.numpy.org/)

### Install

In the main repo directory, you can run python3 setup.py install --user

## Evolve script

In 'project/' you will see scripts names evolve_*.py. These are set up to run simulations in parallel and dump results to 'data/*'

### Recorder class

The recroder class in project/RecordStats.py is a critical piece to this project. It collect's time series data during the simulations. I am currently collecting the following data

* Generation time
* Population size
* Relative load = one minus the mean fitness divided the maximum fitness in the population
* Segregating load = one minus the mean fitness divided the maximum theoretical fitness
* Fixed load = load due to fixations, not considered in other loads
* Number of fixed mutations - neutral and deleterious separately
* Mean number of deleterious mutatiosn per individual
* Sum of derived allele frequencies - neutral and deleterious separately
* Tajima's D at neutral sites
* Tajima's D at all sites
* Heterozygosity (pi) at neutral sites
* Heterozygosity (pi) at all sites
* Fay and Wu's Hprime at neutral sites
* Fay and Wu's Hprime at all sites

## Examples

Please see 'project/Example' for shorter example jupyter notebooks or just take a look at any of the notebooks in the main project directory.
