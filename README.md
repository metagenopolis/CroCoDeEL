# CroCoDeEL : **CRO**ss-sample **CO**ntamination **DE**tection and **E**stimation of its **L**evels

## Installation

> 

## Usage

To run CroCoDeEL, use the following command:

> nextflow run main.nf -w <temporary_file_path> \
--mgsProfilesPath <mgs_profiles> \
--outputDirectoryPath <output_directory_path>

You can add information about the samples plates or information about the samples, particularly if they are longitudinal.

> nextflow run main.nf -w <temporary_file_path> \
--mgsProfilesPath <mgs_profiles_path> \
--outputDirectoryPath <output_directory_path>  \
--plateMapPath  <plate_map_path> \
--longitudinalMetadataPath <longitudinal_metadata_path> 

### 1. Abundance table of MGS
____________________

<center>

|   id_mgs  | sample 1 | sample 2 | sample 3 |    ...   | 
|:----------|:--------:|:--------:|:--------:|:--------:| 
| species 1 |     0    |  7.6e-07 |     0    |    ...   | 
| species 2 |  1.8e-06 | 1.42e-06 |     0    |    ...   | 
|    ...    |    ...   |    ...   |    ...   |    ...   | 

</center>

### 2. Plate map (optional)
____________________

<center>

| samples  |  plate  | 
|:---------|:-------:|
| sample 1 |  plate1 | 
| sample 1 |  plate2 |
| sample 2 |  plate1 |
| sample 3 |  plate2 |
|   ...    |   ...   | 

</center>

### 3. Metadata of longitudinal samples (optional)
____________________

<center>

| samples  |  individual  | 
|:---------|:------------:|
| sample 1 |  individual1 | 
| sample 2 |  individual2 |
| sample 3 |  individual1 |
|   ...    |      ...     | 

</center>

## Citation

## Authors
* Lindsay Goulet: lindsay.goulet@inrae.fr
* Florian Plaza Oñate: florian.plaza-onate@inrae.fr
* Edi Prifti: edi.prifti@ird.fr
* Eugeni Belda: eugeni.belda@ird.fr
* Emmanuelle Le Chatelier: emmanuelle.le-chatelier@inrae.fr
* Guillaume Gautreau: guillaume.gautreau@inrae.fr