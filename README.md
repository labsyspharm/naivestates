# naivestates - Inference of cell states using Naive Bayes

This work is supported by the *NIH Grant 1U54CA225088: Systems Pharmacology of Therapeutic and Adverse Responses to Immune Checkpoint and Small Molecule Drugs* and by the *NCI grant 1U2CCA233262: Pre-cancer atlases of cutaneous and hematologic origin (PATCH Center)*.

**Contributors:** Artem Sokolov

## Running as a Docker container

```
docker run -v /path/to/data/folder:/data labsyspharm/naivestates \
  /app/main.R -i /data/myfile.csv -o /data -m aSMA,CD45,panCK --plots
```

where 

* `-v /path/to/data/folder:/data` maps the local folder containing your data to `/data` inside the container
* `-i /data/myfile.csv` specifies which data file to process
* `-m aSMA,CD45,panCK` specifies the markers of interest (NOTE: comma-delimited, no spaces)
* `--plots` requests the generation of plots showing model fits
