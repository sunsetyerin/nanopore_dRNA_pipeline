Direct RNA seq @ GSC (Oxford Nanopore Technologies)
===================================================

1. [Cloning](#cloning)
2. [Environment](#python-3-environment)
3. [Snakemake usage of this repo](#snakemake-usage-of-this-repo)
4. [Perform the processing of a single sample](#perform-the-processing-of-a-single-sample)
5. [Guppy Basecall](#guppy-basecall)
6. [TailfindR](#tailfindr)
7. [Flair](#flair)
8. [Albacore](#albacore)
9. [EpiNano](#epinano)

## 1. Cloning

Make sure to clone the submodules that comprises dependancies not manageable
by conda (or pip).

```bash
git clone --recursive https://github.com/sunsetyerin/nanopore_dRNA_m6A_pA_tail_isoform.git

or

git clone https://github.com/sunsetyerin/nanopore_dRNA_m6A_pA_tail_isoform.git
cd direct_rna
git submodule update --init --recursive
```

## 2. Python 3 environment

This repository has been developped using python3.7.3 through python3.7.6 in a
[miniconda](https://docs.conda.io/en/latest/miniconda.html) environment.

Though [conda](https://docs.conda.io/projects/conda/en/latest/) is not
required, it is highly recommended in order to manage the dependencies.
Users can choose to manage their dependencies manually and should consult
the `envs/environment.yaml` file as reference.

__Using snakemake to manage the environment__

From a bash terminal:

```bash
snakemake --use-conda --conda-create-envs-only --cores [cores available]
```

__Using conda to install the environment manually.__

```bash
conda env create --file envs/environment.yaml

conda activate direct_rna

snakemake [args]
```

## 3. Snakemake usage of this repo

All scripts, steps and jobs are managed using
[snakemake](https://snakemake.readthedocs.io/en/stable/#)

Input files, workflow parameters and output path are provided through config
files (`configs/*.yaml`).

The workflow is designed to compare direct-RNA samples.

```bash
# To display the jobs and commands including subworkflows
snakemake -np

# to visualize the acyclic rulegraph
snakemake --rulegraph | dot | display -

# To launch the jobs
snakemake --use-conda
```

An example of snakemake profile for "numbers" cluster is available at
`configs/smk_profiles/numbers/config.yaml` and a second one for "slurmgpu" will
be available eventually.

## 4. Perform the processing of a single sample

So far the processing of a sample includes basecalling, polyA tail estimation,
gene/isoform expression and match to Illumina gene expression (if available in
Flair configs, see config file).

Config files contains the values of variables needed to perform analyses.

`configs/sample_config.yaml` is under development and is provided as an 
example. It's meant to be used to apply the analyses on those first data
available from the PromethION platform.

Most of the data directories hosting results will be created *ad-hoc* by the
snakemake manager, although it is possible to create the folders in advance.
**Any of these folders can be replaced by a symlink in order to redirect
voluminous files to appropriate spaces.** The entire `data/` folder can be
written to a project space with a simple:

`ln -s /PATH/TO/SCRATCH/data data`

It is recommended to use symlinks to redirect the data/ folder to a scrath
space since the outputs can get particularly large.

## 5. Guppy Basecall

Oxford Nanopore Technologies apparels uses `MinKNOW` software to basecall
reads in realtime, at the expense of the production of the event table
(formerly produced by `Albacore`). `Albacore` is discontinued and `Guppy`
is currently favored. It can produce a move table during basecall that is
similar to the deprecated event table (make sure to use `--fast5_out`).
`TailfindR` requires one of those tables.

Download available [here](https://community.nanoporetech.com/downloads).

Current configs are using an installation available to dlhost01 & dlhost03.
`/opt/ont-guppy_3.2.2/bin/guppy_basecaller`

## 6. TailfindR

`TailfindR` is designed to be supplied with a single directory containing all
the `.fast5` files. It supports parallelization by specifying the number of
cores to attribute to the `TailfindR` job. However, it was designed and tested
with MinIONs output files and I noticed that it doesn't scale well to a full
PromethION flowcell. The first large usage of `TailfindR` on gphost04 using the
sequencing run of the COLO829 cell line with 1563 `.fast5` files of ~4000 reads
(~6.25M reads) took 1.5 weeks to process all files and failed to produce the
compilation of the results after another 2 weeks.

The chosen workaround to address the issue is to generate temporary directories
each containing a softlink towards a single (or few) `.fast5` file(s) (see
`rule tailfindr_butler`). A `TailfindR` job is then generated for each
temporary directories. Numerous `TailfindR` jobs should be submitted to a
cluster such as `numbers` through `snakemake`.

## 7. Flair

`Flair` is a workflow in itself and requires its own environment to run. Both
the scripts used for direct_rna analysis and environment `.yaml` files are
provided in the `flair` submodule.

Unlike `TailfindR` the input is a single `.fast5` which prompts the need for
catenation of the multiple 4000 reads `.fast5` produced by `Guppy`. 

## 8. Albacore

`Albacore` is a deprecated basecaller that was discontinued by Oxford Nanopore
Technologies to the benefit of the currently maintained `Guppy`. It is included
to this workflow as a requirement for `EpiNano` m6A caller which relies on the
basecall errors of `Albacore` to assign methylated sites.

`EpiNano` states that it requires
[Albacore v2.1.7](https://github.com/enovoa/EpiNano#considerations-when-using-this-software),
to provide accurate prediction. This version is unavailable. The community only
have access to `Albacore v2.3.4`.
[link to Nanoporetech community](https://community.nanoporetech.com/downloads)

*NOTE* `Albacore` is used along with the multiread `.fast5` wrapper available
via the `seamlessf5` package.

## 9. EpiNano

The current version of `EpiNano` allows m6A sites calling at the genomic
position level. The calling at read level is
[under development](https://github.com/enovoa/EpiNano#considerations-when-using-this-software).