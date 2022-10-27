[![DOI](https://zenodo.org/badge/455542613.svg)](https://zenodo.org/badge/latestdoi/455542613)
# P-DOR - a Pipeline to Disentangling Outbreak Rapidly <img src='p-dor_logo.png' align="right" height="159" />

## Introduction
P-DOR is a bioinformatic pipeline for rapid WGS-based bacterial outbreak detection and characterization, carried on by integratic clinical metadata and contextualizing the genomes of interest within a well curated global genomic database.

## Pipeline description

1) The P-DOR framework keeps an updated database of all ESKAPE genomes from the [PATRIC](https://www.patricbrc.org/) collection. Input genomes are joined with a background of the n most similar ones from the database. The selection is performed according to the k-mer distance via Mash.
2) Each genome of the resulting dataset is aligned to a reference genome and the coreSNPs are called. Core-SNPs are defined as single-nucleotide mutated positions flanked by n identical bases in all of the analyzed genomes. Two modes are available:
*Fast mode* is based on the [Mummer](https://github.com/mummer4/mummer) package, aligning each genome to the reference (Nucmer) and SNP detection (show-snp). The coreSNPs calling is performed using an in-house Python script.
 *Extended mode* is based on [ProgressiveMauve](https://darlinglab.org/mauve/user-guide/progressivemauve.html) aligner. The coreSNPs calling is performed using [Purple](https://skynet.unimi.it/index.php/tools/purple-tool/) software. NOTE: Purple is currently under maintainment, please use the fast mode!!!.
3) A Maximum Likelihood phylogeny is inferred using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) software.
4) Epidemiological clusters are assessed on the basis of coreSNPs distances using a threshold value to hypothesize the epidemiological relationship among the strains. This parameter can be set manually (e.g. according to previous studies) or computed by detecting the inflection point in the distribution of the coreSNP distances between all pairs of genomes.
5) A screening for resistance and virulence determinants is also performed through [Abricate](https://github.com/tseemann/abricate).
6) If patient metadata (i.e. ward of hospitalization, date of admission and discharge) are provived, the pipeline reconstruct the route of transmission  through a temporal and spatial representation of the outbreak.

## Installation

```bash
git clone https://github.com/SteMIDIfactory/P-DOR.git
cd P-DOR/
conda env create -f environment.yml
conda activate P-DOR
```
Also, you might need to install P3, the CLI tool to interact with the PATRIC-DB.
You can follow the instructions at https://docs.patricbrc.org/cli_tutorial/cli_installation.html

Quick and dirty code, for Debian users:
```bash
curl -O -L https://github.com/PATRIC3/PATRIC-distribution/releases/download/1.024/patric-cli-1.024.deb
sudo dpkg -i patric-cli-1.024.deb
sudo apt-get -f install
```
A full integration of the P3-API is coming soon

## Quick guide

The command with default settings is:
```bash
python P-DOR.py -q [query genome folder] -sd [Source Dataset sketch file] -ref [reference genome] -snp_thr infl
```
### Options:
```bash
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input data:
  -q <dirname>          query folder containing genomes in .fna format
  -db <filename>         background sketch file - See "Pre-sketched databases" section
  -ref <filename>       reference genome
  -snp_thr SNP_THRESHOLD
                        Threshold number of SNPs to define an epidemic
                        cluster: choices are an integer number or type 'infl'
                        to calculate the threshold by the inflection point of SNPs
                        distribution

Additional arguments:
  -meta META            metadata file; see example file for formatting
                        (default: None)
  -gff [GFF]            annotation file in gff format; if not specified
                        (default) a dummy gff is generated (default: )
  -bkg_folder [BKG_FOLDER]
                        folder containing the genomes from which the
                        background sketch was created; if not specified
                        (default) nearest genomes are downloaded from the
                        PATRIC-DB (default: )
  -call {purple,mummer}
                        Snps calling method (default: mummer)
  -n <int>              Maximum nearest genomes from database (default: 20)
  -t <int>              number of threads (default: 10)

```
### Test run

Test whether the pipeline generates all the expected outputs on your system, run the following command in the P-DOR folder

```
python P-DOR.py -q data/test_query -sd data/sketches.msh -ref data/NJST258.fna -snp_thr infl -n 2 -sd_folder data/test_DB/ -meta data/sample_metadata_table.txt
```

### Pre-sketched databases
Pre-sketched ESKAPE genomes are available at:
```
https://drive.google.com/drive/folders/1lrr0tQn0RRwsHw54zRlZIMIhmdMnZi2Q
```
### Build your own sketch
Run makepdordb.py script indicating the bacterial species. Here, any bacterial species can be indicated.
```
python makepdordb.py -s "Escherichia coli"
```
## Output
1) Summary of resistance and virulence detected.
2) Core-SNPs alignment
3) Core‐SNPs histogram distribution between genome pairs. The dashed bar is set according to the threshold indicating the epidemiological clusters.
4) Heatmap and graph network representing the core-SNPs distances between all pairs of genomes.
5) Maximum Likelihood SNP-based phylogeny with annotated tips according to presence-absence of the genetic determinants of resistance and virulence.   Labels are colored based on the outbreak clusters.
6) Timeline of hospitalized patient and the bacterial samples. The timeline indicates samples isolation based on colonization and infection. The samples are linked according the their core-SNPs distance.

## Coming soon
- Implementation of the [SCOTTI](https://github.com/Taming-the-BEAST/SCOTTI-Tutorial) tool for the reconstruction of the chain of transmission via Bayesian inference.
- Genome assembly, both short and long reads
- A quicker implementation of the Purple algorithm
- Genome characterization: MLST
- Utilization of secondary SNP alignments (e.g. codon 3rd position SNP alignment, intergenic SNP alignment) for downstream analyses

### Citation
DOI: 10.5281/zenodo.6481025
