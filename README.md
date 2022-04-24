# P-DOR <img src='p-dor_logo.png' align="right" height="159" />
Quick and easy outbreak reconstruction pipeline.

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
python P-DOR.py -q [query genome folder] -db [background sketch file] -ref [reference genome] 
```
Options:
```bash
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input data:
  -q <dirname>          query folder containing genomes in .fna format
  -db <dirname>         background sketch file
  -ref <filename>       reference genome

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
  -snp_thr SNP_THR      Threshold number of SNPs to define an epidemic cluster
                        (default: 20)
  -n <int>              Maximum nearest genomes from database (default: 20)
  -t <int>              number of threads (default: 10)

```
### Pre-sketched databases
Pre-sketched ESKAPE genomes are available at:
```
https://drive.google.com/drive/folders/1lrr0tQn0RRwsHw54zRlZIMIhmdMnZi2Q
```
  ### Coming soon
