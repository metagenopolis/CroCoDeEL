# CroCoDeEL : **CRO**ss-sample **CO**ntamination **DE**tection and **E**stimation of its **L**evels

## Introduction

CroCoDeEL is a tool that detects cross-sample (aka well-to-well) contamination in shotgun metagenomic data.\
It accurately identifies contaminated samples but also pinpoints contamination sources and estimates contamination rates.\
CroCoDeEL relies only on species abundance tables and does not need negative controls.

## Installation

CroCoDeEL can be easily installed with conda.

First, clone the git repository in a local directory:
```
git clone https://github.com/metagenopolis/crocodeel.git
cd crocodeel
```

Then, create a specific conda environment with all requirements:
```
conda env create -n crocodeel_env --file conda_env/environment.yml
```

Next, activate the conda environment: 

```
conda activate crocodeel_env
```

Finally, you can test that CroCoDeEL is correctly installed with the following command:
```
python3 crocodeel/crocodeel.py test_install
```

_CroCoDeEL will be released on bioconda soon._

## Quick start
### Input
CroCoDeEL takes as input a species abundance table in TSV format.\
The first column should correspond to species names. The other columns correspond to the abundance of species in each sample.

|   species_name  | sample1 | sample2 | sample3 |    ...   | 
|:----------------|:-------:|:-------:|:-------:|:--------:| 
| species 1       |   0     |  0.05   |   0.07  |    ...   | 
| species 2       |   0.1   |  0.01   |   0     |    ...   | 
|       ...       |   ...   |   ...   |   ...   |    ...   | 

CroCoDeEL works with relative abundances.
The table will automatically be normalized so the abundance of each column equals 1.

**Important**: CroCoDeEL requires the abundance of subdominant species to be accurately estimated.\
We strongly recommend using [the Meteor software suite](https://github.com/metagenopolis/meteor) to generate the species abundance table. 
Otherwise, you can use [sylph](https://github.com/bluenote-1577/sylph).\
We advise against the use of other taxonomic profilers (e.g. MetaPhlan4 or mOTUs) that do not meet this requirement according to our benchmarks.

### Search contamination
Run the following command to search for cross-sample contamination:
```
python3 crocodeel/crocodeel.py search_conta -s species_abundance.tsv -c contamination_events.tsv
```
CroCoDeEL will report all detected contamination events in the _contamination_events.tsv_ output file.\
This TSV file reports for each event the contamination source, the contaminated sample (target) and the estimated contamination rate.\
The score (probability) computed by the Random Forest model as well as species specifically introduced by contamination in the target are also given.

### Visualization of the results
Contaminations events can be visually inspected by generating a PDF file consisting in scatterplots.
```
python3 crocodeel/crocodeel.py plot_conta -s species_abundance.tsv -c contamination_events.tsv -r contamination_events.pdf
```
Each scatterplot compares in a log-scale the species abundance profiles of a contaminated sample (x-axis) and its contamination source (y-axis).\
The contamination line (in red) highlights species specifically introduced by contamination.

### Easy workflow
Alternatively, you can search for cross-sample contamination and create the PDF report in one command.
```
python3 crocodeel/crocodeel.py easy_wf -s species_abundance.tsv -c contamination_events.tsv -r contamination_events.pdf
```

### Results interpretation
CroCoDeEL will probably report false contamination events for samples with similar species abundances profiles (e.g. longitudinal data, animals raised together).\
For non-related samples, CroCoDeEL may occasionally generate false positives that can be filtered out by a human-expert.\
Thus, we strongly recommend inspecting scatterplots of each contamination event to discard potential false positives.\
_We will explain how to do it soon._

## Authors
* Lindsay Goulet: lindsay.goulet@inrae.fr
* Florian Plaza Oñate: florian.plaza-onate@inrae.fr
* Edi Prifti: edi.prifti@ird.fr
* Eugeni Belda: eugeni.belda@ird.fr
* Emmanuelle Le Chatelier: emmanuelle.le-chatelier@inrae.fr
* Guillaume Gautreau: guillaume.gautreau@inrae.fr
