# Author
Ricardo Gomez-Reyes; 2023
# Summary
This markdown summary the process used to curated and enrich sequences from 18S (V9) and taxon-groups related to meiofauna. The steps are as follow:

1. Get database (SILVA 138 SSU Ref NR 99 and pr2_version_4.14.0_SSU)
2. File format convertion (qza inputs for qiime)
3. Dereplicate "unique" sequences (using RESCRIPT in qiime)
4. Make specific amplicon-region dataset  (using RESCRIPT in qiime)
5. File format (Export as tsv and fasta format)
6. Databases processing 
7. Sequences-based (SSU V9) and taxonomic-based (Eukaryota) databases filtering
8. Clean redundant lineage ("_sp", "metagenome", "[a-z]_X", "uncultured", ...)
9. Processed lineage using in-house script (`DB_TAXA2WORMS_CONVERTION*.R`)
10. Filtering relevant Meiofauna groups
11. Merge SILVA and PR2 sequences and worms-lineages files (`CURATED-1389f-1510r_worms_taxonomy.tsv` and `CURATED-1389f-1510r_worms_sequences.fasta`)
12. File format convertion (qza)
13. Dereplication

```bash
conda activate qiime2-2023.2

# 1)
# qiime tools import --show-importable-formats
# qiime tools import --show-importable-formats


# HeaderlessTSVTaxonomyFormat
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

qiime feature-classifier classify-sklearn \
  --i-classifier CLASSIFIER.qza \
  --i-reads dna-sequences.qza \
  --o-classification dna-sequences-taxonomy.qza  

# EXPORT

# Tax
qiime tools export --input-path dna-sequences-taxonomy.qza --output-path CURATED_DB_DIR

```