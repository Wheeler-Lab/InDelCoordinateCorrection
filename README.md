# Snakemake workflow: `InDelCoordinateCorrection`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow to correct genome annotation feature coordinates after polishing a reference genome using short reads from one or more sequencing runs, allowing insertions and deletions. It allows traversing an relationship tree of e.g. mutant cell lines, correcting the genome annotations by traversing down a tree like this:
```
                    Reference
                        |
           ---------------------------
           |                         |
      Cell line 1                Cell line 2
           |                         |
     -------------             -------------
     |           |             |           |
Subclone 1   Subclone 2   Subclone 1   Subclone 2
```
Here, the reference genome is a complete sequence provided as a `.fasta` file together with a corresponding genome annotation `.gff` file. Each of the cell lines lower down in the tree has a set of sequencing read samples that allow polishing of the genome sequence of its parental cell line. The workflow takes care of the polishing and the corresponding shifting of genome annotation feature coordinates. 


## Prerequisites
You need to install `snakemake` to use this workflow. Please refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba). The installation option using `mambaforge` is preferred and well-tested.

You can then go ahead and download the workflow by cloning this repository (`git` needs to be installed to do this):
```bash
git clone https://github.com/Wheeler-Lab/InDelCoordinateCorrection
```

## Usage

In order to use this workflow, customise the given example `config/config.yaml` file for your purposes. You will also need to provide the reference `.gff` and `.fasta` files by putting them into the `resources` folder (`.fasta` files can be automatically downloaded from TriTrypDB if specified correctly). The sequencing sample fastq files need to be provided in the correct directory structure. Please see [here](config/README.md) for more details.

To run the workflow please run
```bash
cd InDelCoordinateCorrection
conda activate snakemake
snakemake --cores all --use-conda --configfile config/config.yaml
```