# Meiofauna
Database curation and data visualization support for Meiofauna metabarcoding data. 

# Author
Ricardo Gomez-Reyes; 2023

# Introduction
Meiofauna is an ecological group defined as an “assemblage of benthic metazoans. It typically the meiofauna of the deep-sea sediments is composed mainly of Nematodes as the most abundant group, reaching proportions many times greater than 90%, followed by harpactidoid copepods and annelids, and by less abundant groups such as Arthropoda, Mollusca, Kinorhyncha, Gastrotricha, Tardígrada, Platyhelminthes, among others (CELIZ, J. A. C. (2018). Meiofaunal biodiversity of deep-sediments from the Gulf of Mexico: a metabarcoding and morphological approach for the establishment of a baseline)

# Summary
This pipeline summarize the process to prepare the database from 18S (V9 subunit) and taxon-groups related to meiofauna. Detailed steps may be consulted in this [doc](https://github.com/RJEGR/Meiofauna/blob/main/RESCRIPT.md). The short steps are as follow:

1. Get database (`SILVA 138 SSU Ref NR 99` `and pr2_version_4.14.0_SSU`)
2. File format convertion (qza inputs for qiime)
3. Dereplicate "unique" sequences (using RESCRIPT in qiime)
4. Make specific amplicon-region dataset  (using RESCRIPT in qiime)
5. File format (Export as tsv and fasta format)
6. Databases processing (in `R`)
7. Sequences-based (SSU V9) and taxonomic-based (Eukaryota) filtering
8. Clean redundant lineage ("_sp", "metagenome", "[a-z]_X", "uncultured", ...)
9. Processed lineage using in-house script (`DB_TAXA2WORMS_CONVERTION*.R`)
10. Filtering relevant Meiofauna groups
11. Merge SILVA and PR2 sequences and worms-lineages files (`CURATED-1389f-1510r_worms_taxonomy.tsv` and `CURATED-1389f-1510r_worms_sequences.fasta`)
12. File format convertion (qza)
13. Dereplication

# Using the curated DB

```bash
conda activate qiime2-2023.2


# 1)
# qiime tools import --show-importable-formats
# qiime tools import --show-importable-formats


# 1)

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format  TSVTaxonomyFormat\
  --input-path CURATED-1389f-1510r_worms_taxonomy.tsv \
  --output-path CURATED-1389f-1510r_worms_taxonomy.qza

# 2)

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-format DNAFASTAFormat \
  --input-path CURATED-1389f-1510r_worms_sequences.fasta \
  --output-path CURATED-1389f-1510r_worms_sequences.qza
  
# 3) Derep 

qiime rescript dereplicate \
    --i-sequences CURATED-1389f-1510r_worms_sequences.qza \
    --i-taxa CURATED-1389f-1510r_worms_taxonomy.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences CURATED-1389f-1510r_worms_derep_sequences.qza \
    --o-dereplicated-taxa CURATED-1389f-1510r_worms_derep_taxonomy.qza
    
# 4) Train class

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads CURATED-1389f-1510r_worms_derep_sequences.qza \
  --i-reference-taxonomy CURATED-1389f-1510r_worms_derep_taxonomy.qza \
  --o-classifier CLASSIFIER.qza

# 5) Classify

WD=/Users/cigom/Documents/MEIOFAUNA_PAPER/INPUTS
cd $WD

qiime feature-classifier classify-sklearn \
  --i-classifier CLASSIFIER.qza \
  --i-reads dna-sequences.qza \
  --o-classification dna-sequences-taxonomy.qza  

# EXPORT

qiime tools export --input-path dna-sequences-taxonomy.qza --output-path CURATED_DB_DIR

# Classi 2

qiime feature-classifier classify-consensus-blast \
  --i-query dna-sequences.qza \
  --i-reference-reads CURATED-1389f-1510r_worms_derep_sequences.qza \
  --i-reference-taxonomy CURATED-1389f-1510r_worms_derep_taxonomy.qza \
  --o-classification dna-sequences-classify-consensus-blast-tax.qza \
  --o-search-results dna-sequences-classify-consensus-blast-results.qza

# EXPORT

qiime tools export --input-path dna-sequences-classify-consensus-blast-tax.qza --output-path classify-consensus-blast_dir


```

# Figures
![Figure 1](https://github.com/RJEGR/Meiofauna/blob/main/Figures/Fig1.png)

![Figure 2](https://github.com/RJEGR/Meiofauna/blob/main/Figures/Fig2.png)
