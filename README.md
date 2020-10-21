# naivestates - Inference of cell states using Naive Bayes

This work is supported by the *NIH Grant 1U54CA225088: Systems Pharmacology of Therapeutic and Adverse Responses to Immune Checkpoint and Small Molecule Drugs* and by the *NCI grant 1U2CCA233262: Pre-cancer atlases of cutaneous and hematologic origin (PATCH Center)*.

# Introduction

`naivestates` is a label-free, cluster-free tool for inferring cell types from quantified marker expression data, based on known marker <-> cell type associations. The tool is designed to be run as a Docker container, but can also be installed in a Conda environment or as an R package. `naivestates` expects as input information about marker expression on a per-cell basis, provided in `.csv` format. One of the columns must contain cell IDs. An example input file may look as follows:

```
CellID,KERATIN,FOXP3,SMA
1,64.18060200668896,193.00334448160535,303.5016722408027
2,54.850202429149796,151.19433198380565,176.3846153846154
3,63.94712643678161,210.43218390804597,483.9448275862069
4,142.01320132013203,227.85808580858085,420.76897689768975
5,56.66379310344828,197.01896551724138,343.7810344827586
6,69.97454545454545,187.59636363636363,267.9709090909091
7,67.57754010695187,185.63368983957218,351.7914438502674
8,64.012,190.02,349.348
9,56.9622641509434,159.79245283018867,236.43867924528303
...
```

# Installation
## Download the container image
Pull the latest version with

```
docker pull labsyspharm/naivestates
```

Alternatively, you can pull a specific version, which is recommended to ensure reproducibility of your analyses. For example, v1.2.0 can be pulled with

```
docker pull labsyspharm/naivestates:1.2.0
```

## Examine the tool usage instructions

```
docker run --rm labsyspharm/naivestates:1.2.0 /app/main.R -h
```

replacing `1.2.0` with the version you are working with. Omit `:1.2.0` entirely if you pulled the latest version above. The flag `--rm` tells Docker to delete the container instance after it finishes displaying the help message.

# Basic usage

At minimum, the tool requires an input file and the list of marker names:

```
docker run --rm -v /path/to/data/folder:/data labsyspharm/naivestates:1.2.0 \
  /app/main.R -i /data/myfile.csv -m aSMA,CD45,panCK
```

where we can make a distinction between Docker-level arguments:

* `--rm` once again cleans up the container instance after it finishes running the code
* `-v /path/to/data/folder:/data` maps the local folder containing your data to `/data` inside the container
* `:1.2.0` specifies the container version that we pulled above

and tool-level arguments:

* `-i /data/myfile.csv` specifies which data file to process
* `-m aSMA,CD45,panCK` specifies the markers of interest (NOTE: comma-delimited, no spaces)

If there is a large number of markers, place their names in a standalone file `markers.txt` with one marker per line. Ensure that the file lives in `/path/to/data/folder/` and modify the Docker call to use the new file:

```
docker run --rm -v /path/to/data/folder:/data labsyspharm/naivestates:1.2.0 \
  /app/main.R -i /data/myfile.csv -m /data/markers.txt
```

## Additional parameters

The following parameters are optional, but may be useful in certain scenarios:

* `-o <path>` - (default: `/data`) Alternative output directory. (Note that any file written to a directory that wasn't mapped with `docker -v` will not persist when the container is destroyed.)
* `--plots <off|pdf|png>` - (default: `off`) Produces QC plots of individual marker fits and summary UMAP plots in .png or .pdf format.
* `--id <name>` - (default: `CellID`) Name of the column that contains cell IDs
* `--log <yes|no|auto>` - (default: `auto`) When a log10 transformation should be applied prior to fitting the data. The tool will do this automatically if it detects large values. Use `--log no` to force the use of original, non-transformed values instead.
* `--sfx <suffix>` - (default: automatically determined) A common suffix on the marker columns (e.g., `_cellMask` or `_nucleiMask`). The suffix will be removed in the output plots and tables to improve readability. Use `$` to force an empty suffix.
* `--umap` - (default: disabled) Include this flag to generate UMAP plots.
* `--mct <filename>` - The tool has a basic marker -> cell type (mct) mapping in `typemap.csv`. More sophisticated mct mappings can be defined by creating a `custom-map.csv` file with two columns: `Marker` and `State`. Ensure that `custom-map.csv` is in `/path/to/data/folder` and point the tool at it with `--mct` (e.g., `/app/main.R -i /data/myfile.csv --mct /data/custom-map.csv -m aSMA,CD45,panCK`)

# Alternative execution environments
## Running in a Conda environment

If you are working in a computational environment that doesn't support Docker, the repository provides a Conda-based alternative. Ensure that `conda` is installed on your system, then 1) clone this repository, 2) instantiate the conda environment and 3) install the tool.

``` bash
git clone https://github.com/labsyspharm/naivestates.git
cd naivestates
conda env create -f conda.yml
conda activate naivestates
R -s -e "devtools::install_github('labsyspharm/naivestates')"
```

The tool can now be used as above by running `main.R`:

``` bash
./main.R -h
./main.R -i /path/to/datafile.csv -m aSMA,CD45,panCK
```

## Running as an R package

The tool can also be installed as an R package directly from GitHub:

``` r
if( !require(devtools) ) install.packages("devtools")
devtools::install_github( "labsyspharm/naivestates" )
```

Example usage:

``` r
library( tidyverse )
library( naivestates )

# Load the original data
X <- read_csv( "datafile.csv" )

# Fit models to channels aSMA, CD45 and panCK
# Specify that cell IDs are in column CellID
GMM <- GMMfit( X, CellID, aSMA, CD45, panCK )

# Plot a fit to one of the markers
plotMarker( GMM, "CD45" )

# Write out the results to results.csv
GMMreshape(GMM) %>% write_csv( "results.csv" )
```
