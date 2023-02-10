# PhageSnake
The PhageSnake was a base bacteriophage(phage) genome analyze protocol, coded by Snakemake.

# Usage
## 0. Installation
#### 0.1 Environment setting
It is easier to rebuild a environment by CONDA.
Used software was list in `./envs/phagesnake.yaml`.

```bash
conda create -n phage_snake -f ./envs/phagesnake.yaml
conda activate phage_snake
conda install snakemake
```

#### 0.2 Database download
Used database would download from [INPHARED](https://github.com/RyanCook94/inphared).
The location of the download database is set in `config.yaml`, and an absolute path is recommended here.
The "db_prefix" was the time of database, e.g. `1Oct2022`.

```bash
snakemake -s ./rules/setup.smk
```

## 1. Input setting
All the input files were phage nucleotide genome assemblies in FASTA type.
Them would be copied into a new folder, default input directory was set as `fna_files`, it could be changed in `config.yaml`.

## 2. Run protocol
When all configs were correct, the protocol could run easily.
```bash
conda activate phage_snake
snakemake -s phagesnake.smk -j 2
```

# Result
## 