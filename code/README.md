# Snakemake workflow: small molecular splicing RNA-seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/20211209_JingxinRNAseq.svg?branch=master)](https://travis-ci.org/snakemake-workflows/20211209_JingxinRNAseq)


## Authors

* Benjamin Fair (@bfairkun)

## Usage

### Step 1: Install workflow and dependencies

If you simply want to use this workflow, clone the [latest release](https://github.com/bfairkun/20211209_JingxinRNAseq).

    git clone git@github.com:bfairkun/20211209_JingxinRNAseq.git

If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

Install snakemake and the workflow's other dependencies via conda/mamba. If conda/mamba isn't already installed, I recommend [installing miniconda](https://docs.conda.io/en/latest/miniconda.html) and then [install mamba](https://github.com/mamba-org/mamba) in your base environment. Then...

    # move to the snakemake's working directory
    cd 20211209_JingxinRNAseq/code
    # Create environment for the snakemake
    mamba env create -f envs/20211209_JingxinRNAseq.yaml
    # And activate the enviroment
    conda activate 20211209_JingxinRNAseq

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`. Use/modify the config yaml files in the `snakemake_profiles/slurm/` profile to run on UChicago RCC Midway with slurm scheduler.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via the included slurm snakemake profile.

    snakemake --profile snakemake_profiles/slurm

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
