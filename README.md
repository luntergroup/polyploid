# polyploid
Benchmarking variant calling in polyploids

### Running experiments

All experiments are reproducable using a Snakemake workflow. First clone the repository:

```shell
$ git clone https://github.com/luntergroup/polyploid && cd polyploid
```

You will need to install Conda, if not already:

```shell
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh # follow instructions, answer 'yes' where asked
$ source ~/.bashrc # assuming you installed conda into your home directory
$ conda update conda
```

Install Snakemake and general dependencies with conda:

```shell
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda create --name polyploid snakemake pysam python-wget
$ conda activate polyploid
```

Each set of experiments is specified in a YAML config file in the `config` directory. They are self-contained, but download links for the [PrecisionFDA Truth v2](https://precision.fda.gov/challenges/10) raw read data must be specified as these are only accessible after authorisation. You can either download the data manually and rename the files appropriatly, or provide the links in a config file, e.g.:

```shell
$ echo "links:\n\tHG002:" >> config/tetraploid_novaseq.yaml
$ echo "\t\t-<link_here>\n\t\t-<link_here>" >> config/tetraploid_novaseq.yaml
```

Then run each experiment as required, e.g:

```shell
$ snakemake --configfile config/tetraploid_novaseq.yaml --use-conda -j 100 --cluster "qsub -cwd -V -P mygroup.prj -q long.qf -pe shmem {threads}"
```


