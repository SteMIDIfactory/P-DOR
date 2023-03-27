[![DOI](https://zenodo.org/badge/455542613.svg)](https://zenodo.org/badge/latestdoi/455542613)
# P-DOR - a Pipeline to Disentangle Outbreaks Rapidly <img src='p-dor_logo.png' align="right" height="159" /> 

## Introduction
P-DOR is a bioinformatic pipeline for rapid WGS-based bacterial outbreak detection and characterization, carried on by integrating clinical metadata and contextualizing the genomes of interest within a well curated global genomic database. 

## Pipeline description

1) The P-DOR framework keeps an updated dataset of all ESKAPE genomes from the [BV-BRC](https://www.bv-brc.org/) collection. This is called Source Dataset (SD) and it is formatted as a single MASH sketch file (.msh). Input genomes are joined with the n most similar ones from the SD to form the Analysis Dataset (AD). The selection is performed according to the k-mer distance via Mash. 
2) Using the [Mummer](https://github.com/mummer4/mummer) package, each genome of the AD is aligned to a reference genome (Nucmer) and mutations are detected (show-snps). Then, an in-house Python script is used to call coreSNPs.
3) A Maximum Likelihood phylogeny is inferred using [iqtree](http://www.iqtree.org/). 
4) Epidemiological clusters are assessed on the basis of coreSNPs distances using a threshold value to hypothesize the epidemiological relationship among the strains. This parameter must be set manually according to previous studies, or estimated according to the distribution of SNP distances among the AD genomes.
5) A screening for resistance and virulence determinants is also performed through [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/).
6) If patient metadata (i.e. ward of hospitalization, date of admission and discharge) are provived, the pipeline reconstruct the route of transmission  through a temporal and spatial representation of the outbreak.

## Installation
###Requirements
- Linux-based OS
- conda 

###Installation command

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
-q is used to indicate the path to the folder containing the genomes that should be analyzed
-sd is the path to the .msh file containing the Source Dataset
-ref is the path to the genome that should be used as reference for alignment. WARNING: in this version of P-DOR the reference genome file must contain only one sequence in fasta format
-snp_thr is the threshold value of the Core-SNP distance between two genomes that belong to the same epidemic cluster

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
Run makepdordb.py script indicating the bacterial species. Here, any bacterial species can be used.
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
Finally, makepdordb can be used to sketch a local collection of genomes (CAUTION: genomes must be in fasta format with the .fna file extension)
  
```
python3.8 makepdordb.py sketch -f [path to the folder containing the local genomes]
```


## Output

The Output files are contained in a folder named "Results_" and date and time when P-DOR was launched. The main Outputs are:

1) Core-SNP alignment (Filename: SNP_alignment.core.fasta)

2) Core‐SNP distribution between genome pairs. The dashed bar corresponds to the threshold set in input (Filename: SNP_frequency_manual_threshold_15_snps.svg, where 15 is the Core-SNP threshold set by the user)

Epidemiological clusters are assessed on the basis of the topology of the phylogeny and of coreSNP distances using a threshold value, which can be set according to the literature or to the SNP-distance distribution in the dataset. Core-SNP distances in bacterial pathogens are known to follow a distribution with multiple peaks. The peak at the lowest SNP number corresponds with the distances within the same epidemic cluster (David et al. 2019). Hence, you should aim to set your SNP threshold value at the inflection point following the first peak.

Here is the Core-SNP distribution after setting the threshold at 15 and running P-DOR for the first time on the dataset.

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/SNP_histogram_first_run.png)
 
In this case the threshold is slightly lower than the actual inflection point.  

After the preliminary run, P-DOR should be run again using a different threshold value, estimated according to the inflection point. In this case, it is around 25 SNPs. 

Here is the Core-SNP distribution histogram, after P-DOR is run the second time

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/SNP_histogram_adjusted.png)


3) Heatmap and graph network representing the core-SNP distances between the query genomes (AD) vs all pairs of genomes (AD+SD) (Filenames: SNP_heatmap_query_vs_all.svg and SNP_clusters_manual_threshold_20_snps.svg)

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/SNP_heatmap_query_vs_all.svg)


4) Maximum Likelihood SNP-based phylogeny with annotated tips according to presence-absence of the genetic determinants of resistance and virulence.   Labels are colored based on the outbreak clusters (Filename: annotated_tree.pdf)

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/annotated_tree_resvir_cluster.svg)

6) Timeline of hospitalized patient and the bacterial samples. The timeline indicates samples isolation based on colonization and infection. The samples are linked according the their core-SNPs distance (Filename: contact_network_plot.svg)

![alt text](https://github.com/SteMIDIfactory/P-DOR/blob/master/output/contact_network_plot.svg)


## Coming soon
- Possibility to use a reference genome that is split into multiple contigs/chromosomes
- Implementation of the [SCOTTI](https://github.com/Taming-the-BEAST/SCOTTI-Tutorial) tool for the reconstruction of the chain of transmission via Bayesian inference
- Genome assembly, both short and long reads
- Genome characterization: MLST
- Utilization of secondary SNP alignments (e.g. codon 3rd position SNP alignment, intergenic SNP alignment) for downstream analyses

### Citation
DOI: 10.5281/zenodo.6481025
