[![DOI](https://zenodo.org/badge/455542613.svg)](https://zenodo.org/badge/latestdoi/455542613)
# P-DOR - a Pipeline to Disentangle Outbreaks Rapidly <img src='p-dor_logo.png' align="right" height="159" /> 

## Introduction
P-DOR is a bioinformatic pipeline for rapid WGS-based bacterial outbreak detection and characterization, carried on by integrating clinical metadata and contextualizing the genomes of interest within a well curated global genomic database. 

## Pipeline description

1) The P-DOR framework keeps an updated dataset of all ESKAPE genomes from the [BV-BRC](https://www.bv-brc.org/) collection. This is called Source Dataset (SD). Input genomes are joined with the n most similar ones from the SD to form the Analysis Dataset (AD). The selection is performed according to the k-mer distance via Mash. 
2) Using the [Mummer](https://github.com/mummer4/mummer) package, each genome of the AD is aligned to a reference genome (Nucmer) and mutations are detected (show-snps). Then, an in-house Python script is used to call coreSNPs.
3) A Maximum Likelihood phylogeny is inferred using [iqtree](http://www.iqtree.org/). 
4) Epidemiological clusters are assessed on the basis of coreSNPs distances using a threshold value to hypothesize the epidemiological relationship among the strains. This parameter must be set manually according to previous studies, or estimated according to the distribution of SNP distances among the AD genomes.
5) A screening for resistance and virulence determinants is also performed through [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/).
6) If patient metadata (i.e. ward of hospitalization, date of admission and discharge) are provived, the pipeline reconstruct the route of transmission  through a temporal and spatial representation of the outbreak.

## Installation

```bash
git clone https://github.com/SteMIDIfactory/P-DOR.git
cd P-DOR/
conda env create -f environment.yml
conda activate P-DOR
```


## Quick guide

The command with default settings is:
```bash
python3.8 P-DOR.py -q [query genome folder] -sd [background sketch file] -ref [reference genome] -snp_thr 20 
```
### Options:
```
optional arguments:
  -h, --help            show this help message and exit
  -TEST                 TEST MODE. Warning: this option overrides all other arguments (default: False)
  -v, --version         show program's version number and exit

Input data:
  -q <dirname>          query folder containing genomes in .fna format (default: None)
  -sd <dirname>         Source Dataset (SD) sketch file (default: None)
  -ref <filename>       reference genome (default: None)
  -snp_thr <int>        Threshold number of SNPs to define an epidemic cluster (default: None)

Additional arguments:
  -amrf                 if selected, P-DOR uses amrfinder-plus to search for antimicrobial resistance and virulence genes in the entire Analysis Dataset (default: False)
  -meta META            metadata file; see example file for formatting (default: None)
  -sd_folder [BKG_FOLDER]
                        folder containing the genomes from which the Source Dataset (SD) sketch was created (default: )
  -borders <int>        length of the regions at the contig extremities from which SNPs are not called (default: 20)
  -min_contig_length <int>
                        minimum contig length to be retained for SNPs analysis (default: 500)
  -snp_spacing <int>    Number of bases surrounding the mutated position to call a SNP (default: 10)
  -species SPECIES      full species name of the genomes in the query folder. Needs to be written within single quotes, e.g.: -species 'Klebsiella pneumoniae'. This option is used for a quicker and more precise gene search with AMRFinderPlus (default: )
  -n <int>              Maximum closest genomes from Source Dataset (SD) (default: 20)
  -t <int>              number of threads (default: 2)
```
### Test run

Test whether the pipeline generates all the expected outputs on your system, run the following command. 

```
python3.8 P-DOR.py -TEST
```

### Pre-sketched databases
Pre-sketched ESKAPE genomes are available at:
```
https://drive.google.com/drive/folders/1lrr0tQn0RRwsHw54zRlZIMIhmdMnZi2Q
```
### Build your own sketch
Run makepdordb.py script indicating the bacterial species. Here, any bacterial species can be indicated.
```
python3.8 makepdordb.py download -s "Escherichia coli" 
```
You can also use makepdordb.py to update a pre-existing collection
```
python3.8 makepdordb.py download -s "Escherichia coli" -f [path to the pre-existing folder]
```
makepdordb.py uses assembly-stats to check which genomes of the BV-BRC database are already in your folder. This process can be performed using multiple threads (default= 1)

```
python3.8 makepdordb.py download -s "Escherichia coli" -f [path to the pre-existing folder] -t [number of threads to use when checking genomes]
```
Finally, makepdordb can be used to sketch a local collection of genomes
  
```
python3.8 makepdordb.py sketch -f [path to the folder containing the local genomes]
```


## Output
1) Summary of resistance and virulence detected.

2) Core-SNPs alignment

3) Core‐SNPs histogram distribution between genome pairs. The dashed bar is set according to the threshold indicating the epidemiological clusters.

Epidemiological clusters are assessed on the basis of the topology of the phylogeny and of coreSNP distances using a threshold value, which can be set according to the literature or to the SNP-distance distribution in the dataset (inflection point). 

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/SNP_histogram_first_run.png)
 
In this case the threshold was set to 15 SNPs, which is slightly less than the actual inflection point.  

Here, the threshold can be adjusted after the preliminary run of the pipeline according to the inflection point that is around 25 SNPs, a more indicated threshold to be set in the analysis. 

Once the adjusted threshold is estimated, P-DOR can be run again to obtain a better assessment of the epidemiological clusters.

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/SNP_histogram_first_run.png)


4) Heatmap and graph network representing the core-SNPs distances between all pairs of genomes.

5) Maximum Likelihood SNP-based phylogeny with annotated tips according to presence-absence of the genetic determinants of resistance and virulence.   Labels are colored based on the outbreak clusters.

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/annotated_tree_resvir_cluster.svg)


6) Timeline of hospitalized patient and the bacterial samples. The timeline indicates samples isolation based on colonization and infection. The samples are linked according the their core-SNPs distance.

## Coming soon
- Implementation of the [SCOTTI](https://github.com/Taming-the-BEAST/SCOTTI-Tutorial) tool for the reconstruction of the chain of transmission via Bayesian inference.
- Genome assembly, both short and long reads
- Genome characterization: MLST
- Utilization of secondary SNP alignments (e.g. codon 3rd position SNP alignment, intergenic SNP alignment) for downstream analyses

### Citation
DOI: 10.5281/zenodo.6481025
