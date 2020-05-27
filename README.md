# naivestates - Inference of cell states using Naive Bayes

This work is supported by the *NIH Grant 1U54CA225088: Systems Pharmacology of Therapeutic and Adverse Responses to Immune Checkpoint and Small Molecule Drugs* and by the *NCI grant 1U2CCA233262: Pre-cancer atlases of cutaneous and hematologic origin (PATCH Center)*.

**Contributors:** Artem Sokolov

# Running as a Docker container

## Download the container image
Pull the latest version with

```
docker pull labsyspharm/naivestates
```

Alternatively, you can pull a specific version, which is recommended to ensure reproducibility of your analyses. For example, v1.1.0 can be pulled with

```
docker pull labsyspharm/naivestates:1.1.0
```

## Examine the tool usage instructions

```
docker run --rm labsyspharm/naivestates:1.1.0 /app/main.R -h
```

replacing `1.1.0` with the version you are working with. Omit `:1.1.0` entirely if you pulled the latest version above. The flag `--rm` tells Docker to delete the container instance after it finishes displaying the help message.

## Apply the tool to quantification data that captures cell-level marker expression

```
docker run --rm -v /path/to/data/folder:/data labsyspharm/naivestates:1.1.0 \
  /app/main.R -i /data/myfile.csv -o /data -m aSMA,CD45,panCK --plots --log auto --id CellID
```

where we can make a distinction between Docker-level arguments:

* `--rm` once again cleans up the container instance after it finishes running the code
* `-v /path/to/data/folder:/data` maps the local folder containing your data to `/data` inside the container
* `:1.1.0` specifies the container version that we pulled above

and tool-level arguments:

* `-i /data/myfile.csv` specifies which data file to process
* `-o /data` tells the tool to write output to `/data` inside the container (recall that this is mapped to `/path/to/data/folder` above)
* `-m aSMA,CD45,panCK` specifies the markers of interest (NOTE: comma-delimited, no spaces)
* `--plots` requests the generation of plots showing model fits
* `--log` can be one of `<yes|no|auto>` (where `auto` is the default), which specifies whether the tool should apply a log10 transformation prior to fitting the data
* `--id` tells the tool which column contains Cell IDs. If omitted, the tool will look for `CellID`.

### Fitting a large number of markers

If there is a large number of markers, place their names in a standalone file `markers.txt` with one marker per line. Ensure that the file lives in `/path/to/data/folder/` and modify the above Docker call to use the new file:

```
docker run --rm -v /path/to/data/folder:/data labsyspharm/naivestates:1.1.0 \
  /app/main.R -i /data/myfile.csv -o /data -m /data/markers.txt
```

Note that the optional parameters `--plot`, `--log` and `--id` were omitted for clarity.
