#####
# Metabarcoding pipeline for
# Shark-dust: High-throughput DNA sequencing of processing residues unveils rife trade in endangered sharks and rays
#
# Andhika P. Prasetyo,1,2* Joanna M. Murray,3 Firdaus Agung,4 Efin Muttaqin,5 Naiara G. Sales,6 Allan D. McDevitt,1,7† Stefano Mariani8†
#
# *Corresponding author. Email: a.p.prasetyo@edu.salford.ac.uk
#####

# Working directory
cd sharkdust/Data/0_raw_data

### 1. Paired-end alignment: merges R1 and R2, then splits in "Good" and "Bad" based on an alignment score cut-off of 40.

illuminapairedend -r Dust_S1_L001_R1_001.fastq Dust_S1_L001_R2_001.fastq | obiannotate -S goodali:'"Good_Dust" if score>40.00 else "Bad_Dust"' | obisplit -t goodali

## Output:
# Good_Dust.fastq (9,366,531 reads)
# Bad_Dust.fastq (185,682 reads)  


### 2. Demultiplex (need ngsfilter formatted .tsv file)
ngsfilter -t 1_ngsfilter_Dust.txt --fasta-output -u unidentified_Dust.fasta Good_Dust.fastq > Dust.filtered.fasta

## Output:
# Dust.filtered.fasta (8,075,259 reads)
# unidentified_Dust.fasta (3,091,480 reads)   

### 3. Filter the seqs with length between 140 and 190 bp and remove seqs with 'N'
obigrep -p 'seq_length>140' -p 'seq_length<190' -s '^[ACGT]+$' Dust.filtered.fasta > Dust.filtered_length.fasta
 
## Output:
# Dust.filtered_length.fasta (7,990,012 reads)

### 4. Group the unique seqs (obiuniq)
obiuniq -m sample Dust.filtered_length.fasta > Dust.unique.fasta

## Output:
# Dust.unique.fasta (718,342 reads)

### 5. Exchange the identifier to a short index _(obiannotate)_
obiannotate --seq-rank Dust.unique.fasta | obiannotate --set-identifier '"'Dust'_%09d" % seq_rank' > Dust.fasta

## Output:
# Dust.fasta (718,342 reads)

### 6. Change formats to vsearch (in R_scripts_metabarpark-master)
Rscript /home/naiara/SOFTWARE/scripts/owi_obifasta2vsearch -i Dust.fasta -o Dust.vsearch.fasta

## Output:
# Dust.vsearch.fasta (587,398 reads)

### 7. Searching for chimeras
### 7.1. Edit input file before run vsearch. By removing space in every header by changing “ ;size=” by “;size=”. 
edit Dust.vsearch.fasta
vim Dust.vsearch.fasta
:%s/ ;size=/;size=/g
:wq

### 7.2. Run UCHIME de novo in VSEARCH 
vsearch --uchime_denovo Dust.vsearch.fasta --sizeout --minh 0.90 --nonchimeras Dust.nonchimeras.fasta --chimeras Dust.chimeras.fasta --uchimeout Dust.uchimeout.txt

## Output:
# vsearch v2.13.3_linux_x86_64, 503.8GB RAM, 24 cores
# https://github.com/torognes/vsearch

# Reading file Dust.vsearch.fasta 100%
# 106396997 nt in 587738 seqs, min 141, max 189, avg 181
# Masking 100%
# Sorting by abundance 100%
# Counting k-mers 100%
# Detecting chimeras 100%
# Found 9676 (1.6%) chimeras, 578062 (98.4%) non-chimeras,
# and 0 (0.0%) borderline sequences in 587738 unique sequences.
# Taking abundance information into account, this corresponds to
# 15232 (0.2%) chimeras, 6457063 (99.8%) non-chimeras,
# and 0 (0.0%) borderline sequences in 6472295 total sequences.

### 8. Clustering using SWARM
swarm -d 1 -z -t 10 -o Dust_d1_SWARM3nc_output -s Dust_d1_SWARM3nc_stats -w Dust_d1_SWARM3nc_seeds.fasta Dust.nonchimeras.fasta

## Output:
# CPU features:      mmx sse sse2 sse3 ssse3 sse4.1 sse4.2 popcnt avx avx2
# Database file:     Dust.nonchimeras.fasta
# Output file:       Dust_d1_SWARM3nc_output
# Statistics file:   Dust_d1_SWARM3nc_stats
# Resolution (d):    1
# Threads:           10
# Break OTUs:        Yes
# Fastidious:        No
# 
# Reading database:  100%  
# Indexing database: 100%  
# Database info:     106396997 nt in 587738 sequences, longest 189 nt
# Hashing sequences: 100%  
# Clustering:        100%  
# Writing swarms:    100%  
# Writing seeds:     100%  
# Writing stats:     100%  
# 
# Number of swarms:  93158
# Largest swarm:     331443
# Max generations:   47

### 9. Recount after SWARM (in R_scripts_metabarpark-master)
owi_recount_swarm Dust_d1_SWARM3nc_output Dust.new.tab

## Output:
# Reading swarm database...
# Cluster database read including 93158 total clusters.
# Calculating number of reads in each cluster
# Kept only 1733 clusters of size greater than or equal to 2 reads.
# Reading tabulated database. This could take a while...
# Database read including 587738 total different sequences and 40 samples.
# Kept only 496313 sequences for calculations.
# File Dust_d1_SWARM3nc_output.counts.csv writtens...

### 10 Select only non singleton MOTUs
### 10.1. Edit input file before Select only non singletons. By adding a space in every header by changing “;size=” by “; size=”. 
vim  Dust_d1_SWARM3nc_seeds.fasta
:%s/;size=/; size=/g
:wq

## Output:
# Dust_d1_SWARM3nc_seeds.fasta (46,579 reads)

### 10.2. Select only non singletons
obigrep -p 'size>1' Dust_d1FORMAT_SWARM3nc_seeds.fasta > Dust_d1FORMAT_seeds_nonsingletons.fasta

## Output:
# Dust_d1_seeds_nonsingletons.fasta (2,065 reads)

### 10.3. Annotate with ecotag
ecotag -d Taxo_May2018/Eukarya -R db_Eukarya_Miya_2020.fasta Dust_d1FORMAT_seeds_nonsingletons.fasta > Dust_d1FORMAT_ecotag_db2020.fasta 

## Output:
# 82.698% of the alignments was cached
# Dust_d1FORMAT_ecotag_db2020.fasta  (2,065 reads)

### 11. Add high level taxa (in R_scripts_metabarpark-master)
Rscript /home/andhika/extrastorage/eDNA_scripts/owi_add_taxonomy Dust_d1_ecotag_db2020.fasta

## Output:
# Reading ecotagged fasta file
# Read 1733 records
# Output file Dust_d1_ecotag_db2020.fasta.annotated.csv written with 1733 sequences

### 12. Combine Ecotag and abundance files (in R_scripts_metabarpark-master)
Rscript owi_combine -i Dust_d1_ecotag_db2020.fasta.annotated.csv -a Dust_d1_SWARM3nc_output.counts.csv -o 2_Dust_d1_MOTUs_db2020.csv

## Output:
# Reading ecotag database...
# Ecotag database read including 1733 total MOTUs.
# Reading abundance database...
# Abundances database read including 1733 total MOTUs and 40 samples.
# File Dust_d1_MOTUs_db202.csv written, including 1733 MOTUs with 6,380,200 total reads in 40 samples.
# (1733 non-singletons MOTUs).
