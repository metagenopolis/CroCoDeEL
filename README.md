#  CroCoDeEL : **CRO**ss-sample **CO**ntamination **DE**tection and **E**stimation of its **L**evel 🐊

[![install with conda](https://img.shields.io/conda/vn/bioconda/crocodeel?color=green&label=bioconda%2Fcrocodeel&logo=anaconda)](https://anaconda.org/bioconda/crocodeel)
[![PyPI](https://img.shields.io/pypi/v/crocodeel?label=pypi%20package)](https://pypi.org/project/crocodeel/)

## Introduction

CroCoDeEL is a tool that detects cross-sample contamination (aka well-to-well leakage) in shotgun metagenomic data.\
It accurately identifies contaminated samples but also pinpoints contamination sources and estimates contamination rates.\
CroCoDeEL relies only on species abundance tables and does not need negative controls nor sample position during processing (i.e. plate maps).

<p align="center">
  <img src="docs/logos/logo.webp" width="350" height="350" alt="logo">
</p>

## Installation

CroCoDeEL is available on bioconda:
```
conda create --name crocodeel_env -c conda-forge -c bioconda crocodeel
conda activate crocodeel_env
```

Alternatively, you can use pip:
```
pip install crocodeel
```

Docker and Singularity containers are also available on [BioContainers](https://biocontainers.pro/tools/crocodeel)

## Installation test

You can test that CroCoDeEL is correctly installed with the following command:
```
crocodeel test_install
```

## Quick start
### Input
CroCoDeEL takes as input a species abundance table in TSV format.\
The first column should correspond to species names. The other columns correspond to the abundance of species in each sample.
An example is available [here](crocodeel/test_data/mgs_profiles_test.tsv).

|   species_name  | sample1 | sample2 | sample3 |    ...   | 
|:----------------|:-------:|:-------:|:-------:|:--------:| 
| species 1       |   0     |  0.05   |   0.07  |    ...   | 
| species 2       |   0.1   |  0.01   |   0     |    ...   | 
|       ...       |   ...   |   ...   |   ...   |    ...   | 

CroCoDeEL works with relative abundances.
The table will automatically be normalized so the abundance of each column equals 1.

**Important**: CroCoDeEL requires the abundance of subdominant species to be accurately estimated.\
We strongly recommend using [the Meteor software suite](https://github.com/metagenopolis/meteor) to generate the species abundance table.\
Alternatively, you can use [sylph](https://github.com/bluenote-1577/sylph) although low-level contaminations may go unnoticed.\
We advise against the use of other taxonomic profilers (e.g. MetaPhlan4 or mOTUs) that do not meet this requirement according to our benchmarks.

### Search contamination
Run the following command to search for cross-sample contamination:
```
crocodeel search_conta -s species_abundance.tsv -c contamination_events.tsv
```
CroCoDeEL will report all detected contamination events in the _contamination_events.tsv_ output file.\
This TSV file reports for each event the contamination source, the contaminated sample (target) and the estimated contamination rate.\
The score (probability) computed by the Random Forest model as well as species specifically introduced by contamination in the target are also given.
An example is available [here](crocodeel/test_data/results/contamination_events.tsv).

### Visualization of the results
Contaminations events can be visually inspected by generating a PDF file consisting in scatterplots.
```
crocodeel plot_conta -s species_abundance.tsv -c contamination_events.tsv -r contamination_events.pdf
```
Each scatterplot compares in a log-scale the species abundance profiles of a contaminated sample (x-axis) and its contamination source (y-axis).\
The contamination line (in red) highlights species specifically introduced by contamination.
An example is available [here](crocodeel/test_data/results/contamination_events.pdf).

### Easy workflow
Alternatively, you can search for cross-sample contamination and create the PDF report in one command.
```
crocodeel easy_wf -s species_abundance.tsv -c contamination_events.tsv -r contamination_events.pdf
```

### Results interpretation
CroCoDeEL will probably report false contamination events for samples with similar species abundances profiles (e.g. longitudinal data, animals raised together).\
For non-related samples, CroCoDeEL may occasionally generate false positives that can be filtered out by a human-expert.\
Thus, we strongly recommend inspecting scatterplots of each contamination event to discard potential false positives.\
Please check the [wiki](https://github.com/metagenopolis/CroCoDeEL/wiki) for more information.
