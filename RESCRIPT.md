# Rescript
...
Installation
```bash

# 1) Install or upgrade QIIME (OS X Intel and/or Linux)
conda update conda

conda install wget

wget https://data.qiime2.org/distro/core/qiime2-2023.2-py38-osx-conda.yml

conda env create -n qiime2-2023.2 --file qiime2-2023.2-py38-osx-conda.yml

rm qiime2-2023.2-py38-osx-conda.yml

conda activate qiime2-2023.2

qiime --help

# 2) Install rescript pluging or/and dependencies

conda activate qiime2-2023.2

pip install git+https://github.com/bokulich-lab/RESCRIPt.git

# Now rescript pluging will be activaded (Pipeline for reference sequence annotation and curation)
# To activate this environment, use
#
#     $ conda activate qiime2-2023.2
#
# To deactivate an active environment, use
#
#     $ conda deactivate

```

## Download DB
https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494

### 1) SILVA
The SILVA 138 SSU Ref NR 99 database has been downloaded using rescript plugging for subsequent data analysis (*SSURef_Nr99_tax_silva.fasta),

```bash
WD=/Users/cigom/Documents/MEIOFAUNA_PAPER/RESCRIPT

cd $WD

conda activate qiime2-2023.2

qiime rescript --help

# 1) Download db and target rannks: domain phylum class order family genus

target="SSURef_NR99"
version="138.1"

qiime rescript get-silva-data \
    --p-version $version \
    --p-target $target \
    --p-no-rank-propagation \
    --p-ranks domain phylum class order family genus \
    --p-include-species-labels \
    --o-silva-sequences ${target}-${version}-rna-seqs.qza \
    --o-silva-taxonomy ${target}-${version}-tax.qza
    
    
#Saved FeatureData[RNASequence] to: SSURef_NR99-138.1-rna-seqs.qza
#Saved FeatureData[Taxonomy] to: SSURef_NR99-138.1-tax.qza
```
### NCBI
https://forum.qiime2.org/t/using-rescript-to-compile-sequence-databases-and-taxonomy-classifiers-from-ncbi-genbank/15947
- https://www.ncbi.nlm.nih.gov/refseq/targetedloci/

SELECTED GENBANK sequences (Animals, Plants and Protists)
NCBI makes separate databases available for download. In this case, "nr" is non-redundant protein, "nt" is **non-redundant nucleotide**

For  unknown-ASV classification sequences may be screened against the NCBI nr/nt and databases using DIAMOND BLASTX . Further verification of the quality of the rRNA sequences may be performed using the Ribovore

```bash
# qiime rescript get-ncbi-data --help

mkdir -p NCBI_nt_db
cd NCBI_nt_db
echo `date +%Y-%m-%d` > download_date.txt

# wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*" # ~ 20 chunks from db

# 18S ribosomal RNA sequences (SSU) from Fungi type and reference material (BioProject PRJNA39195) [blast db files <makeblast>]

18S_fungal_sequences.tar.gz # It is the V4 and V5

# Small subunit ribosomal RNA sequences for eukaryotic sequences
SSU_eukaryote_rRNA.tar.gz
# Additional taxonomy information for the databases 

taxdb.tar.gz           
1758664

# non-redundant protein sequence database with entries from GenPept, Swissprot, PIR, PDF, PDB, and RefSeq

nr.gz*

# nucleotide sequence database, with entries from all traditional divisions of GenBank, EMBL, and DDBJ; excluding bulk divisions (gss, sts, pat, est, htg)

nt.gz*  

# Downloading 24/03/23

# wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz" # 252 Gb ?

targz_files=$(find . | grep \.tar\.gz$ | sed 's/\.\///g')
for f in $targz_files; do tar -xzvf $f; done
rm $targz_files

# Como actualizar NCBI db frequently https://chk.ipmb.sinica.edu.tw/wiki/doku.php/computers/bioinfo_server_configuration#ncbi_db

```
Then blast
```bash
# 7. Formatting a FASTA file into a BLASTable database

# FASTA files need to be formatted with makeblastdb before they can be used in local blast search. For those from NCBI, the following makeblastdb commands are recommended:

For nucleotide fasta file:   makeblastdb -in input_db -dbtype nucl -parse_seqids
For protein fasta file:      makeblastdb -in input_db -dbtype prot -parse_seqids

In general, if the database is available as BLAST database, it is better to use the 
preformatted database.
```

### Ribotyper
https://github.com/ncbi/ribovore/blob/master/documentation/ribotyper.md#top
```bash
ribotyper -f $RIBOSCRIPTSDIR/testfiles/example-16.fa test
```
### PR2
Ranks:
domain	rank 1
supergroup	rank 2
division	rank 3
class	rank 4
order	rank 5
family	rank 6
genus	rank 7
species	Assigned species - rank 8

```bash
# 1) Sequences in fasta format with accession in description line
wget https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_mothur.fasta.gz

# 2) Taxonomy of each sequence separated from the accession number by a tabulation
wget https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_mothur.tax.gz

# 3)

gunzip *gz

#4) Import

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path pr2_version_4.14.0_SSU_mothur.tax \
  --output-path pr2_version_4.14.0_SSU_tax.qza

# 4.1)

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-format DNAFASTAFormat \
  --input-path pr2_version_4.14.0_SSU_mothur.fasta \
  --output-path pr2_version_4.14.0_SSU_seqs.qza
  
  
qiime rescript evaluate-seqs \
    --i-sequences pr2_version_4.14.0_SSU_seqs.qza \
    --o-visualization pr2_version_4.14.0_SSU_seqs.qzv
```
#### PR2: Prepare specific amplicon region DB
```bash
# 1)
qiime rescript dereplicate \
    --i-sequences pr2_version_4.14.0_SSU_seqs.qza \
    --i-taxa pr2_version_4.14.0_SSU_tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences pr2_version_4.14.0_SSU-seqs-euk-derep-uniq.qza \
    --o-dereplicated-taxa pr2_version_4.14.0_SSU_tax-derep-uniq.qza
    
# 2)

qiime feature-classifier extract-reads \
    --i-sequences pr2_version_4.14.0_SSU-seqs-euk-derep-uniq.qza \
    --p-f-primer TTGTACACACCGCCC \
    --p-r-primer CCTTCYGCAGGTTCACCTAC \
    --p-identity 0.7 \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads pr2_version_4.14.0_SSU-seqs-euk-derep-uniq-1389f-1510r.qza
    
# Saved FeatureData[Sequence] to: pr2_version_4.14.0_SSU-seqs-euk-derep-uniq.qza
# Saved FeatureData[Taxonomy] to: pr2_version_4.14.0_SSU_tax-derep-uniq.qza

# IMPORT
qiime tools export --input-path pr2_version_4.14.0_SSU-seqs-euk-derep-uniq-1389f-1510r.qza --output-path pr2_version_4.14.0_SSU

# 2) Fasta
qiime tools export --input-path SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.qza --output-path SSURef_NR99-138.1

    
```

## (Optional) Evaluate seqs and tax
```bash
# 1)
qiime rescript evaluate-seqs \
    --i-sequences SSURef_NR99-138.1-dna-seqs.qza \
    --o-visualization SSURef_NR99-138.1-dna-seqs.qzv
    
# 2)
qiime rescript evaluate-taxonomy \
    --i-taxonomies SSURef_NR99-138.1-tax.qza \
    --o-taxonomy-stats SSURef_NR99-138.1-tax.qzv

```

## 2) FeatureData conversion (rna to dna)
If you'd like to be able to jump to steps that only take FeatureData[Sequence] as input you can convert your data to FeatureData[Sequence] like so:

```bash
qiime rescript reverse-transcribe \
    --i-rna-sequences SSURef_NR99-138.1-rna-seqs.qza \
    --o-dna-sequences SSURef_NR99-138.1-dna-seqs.qza
```

## Optional cleaning
Es algo de truco esta parte, nos interesa primero flankear el fragmento V9 (Paso "Make specific amplicon ...") y despues evaluar secuencias con nucleotidos redundantes (`cull-seqs`), remover secuencias cortas (`filter-seqs-length-by-taxon`) y grupos taxonomicos relevantes.

###  (Optional) Filter out low-quality sequences
```bash
qiime rescript cull-seqs \
    --i-sequences SSURef_NR99-138.1-dna-seqs.qza \
    --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza
```

### Filter DB by Domain (not working properly)

```bash
# Filtering sequences by taxonomy
# One or more of the following filter settings are required: min_lens, max_lens.

# Using --p-min-lens to allow V9 fragments

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences SSURef_NR99-138.1-dna-seqs.qza \
    --i-taxonomy SSURef_NR99-138.1-tax.qza \
    --p-labels Eukaryota \
    --p-min-lens 900 \
    --o-filtered-seqs SSURef_NR99-138.1-dna-seqs-euk.qza \
    --o-discarded-seqs SSURef_NR99-138.1-dna-discarded-seqs-euk.qza
    
# 2) 

qiime rescript filter-taxa \
    --i-taxonomy SSURef_NR99-138.1-tax.qza \
    --p-include Eukaryota \
    --o-filtered-taxonomy SSURef_NR99-138.1-tax-euk.qza
```

## 3) Dereplicate sequences
Using the default `uniq` approach. It will retain identical sequence records that have differing taxonomies ...

```bash
qiime rescript dereplicate \
    --i-sequences SSURef_NR99-138.1-dna-seqs.qza \
    --i-taxa SSURef_NR99-138.1-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences SSURef_NR99-138.1-dna-seqs-derep-uniq.qza \
    --o-dereplicated-taxa SSURef_NR99-138.1-tax-derep-uniq.qza
    
# Saved FeatureData[Sequence] to: SSURef_NR99-138.1-dna-seqs-euk-derep-uniq.qza
# Saved FeatureData[Taxonomy] to: SSURef_NR99-138.1-tax-derep-uniq.qza
```

## 4) Make specific amplicon-region  DB
generate an amplicon-specific DB  for more robust taxonomic classification of your data (Werner et al. 2011, Bokulich et al. 2018).

(From the manuscript) Following the Earth Microbiome Project, the V9 region (~260 pb) of the minor subunit of the ribosomal RNA gene (18S SSU rRNA) was amplified with Euk_1389f (TTGTACACACCGCCC) and Euk_1510r (CCTTCYGCAGGTTCACCTAC) primers. Sequencing of amplicons was completed on an Illumina HiSeq 2500 instrument to produce 250 bp paired-end reads. Roche adapters (shown in bold) and the “barcodes” or 5-base keys (shown as X's) for distinguishing between samples on a single 454 run are detailed in the table. Our 1380F, 1389F and 1510R V9 primers were synthesized at Invitrogen (Carlsbad, CA), HPLC or cartridge purified and engineered with 5 base keys, avoiding the use of a C at the terminal position of the key preceding the 1380F primer.


```bash
# The universal-specific forward/eukaryotic-specific reverse primer combination: 1389F/1510R (Amaral-Zettler et al 2009)

#- Euk_1389f (Universal):  GCCTCCCTCGCGCCATCAGXXXXXTTGTACACACCGCCC
#- Euk_1510r: GCCTTGCCAGCCCGCTCAGCCTTCYGCAGGTTCACCTAC

# the primer pair are botu, 1389 F 5′-TTGTACACACCGCCC-3′ and 1510 R 5′-CCTTCYGCAGGTTCACCTAC-3′ 
# 
qiime feature-classifier extract-reads \
    --i-sequences SSURef_NR99-138.1-dna-seqs-derep-uniq.qza \
    --p-f-primer TTGTACACACCGCCC \
    --p-r-primer CCTTCYGCAGGTTCACCTAC \
    --p-identity 0.7 \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r.qza
    
# 2)


qiime rescript evaluate-seqs \
    --i-sequences SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.qza \
    --o-visualization SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.qzv
    
# IMPORT V9 FILE

# FASTA

qiime tools export --input-path SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r.qza --output-path SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r_dir

# 2) TAXA
qiime tools export --input-path SSURef_NR99-138.1-tax-derep-uniq.qza --output-path SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r_dir

```

## 5) Train the Classify (sklearn)
```bash
ln -s /Users/cigom/Documents/MEIOFAUNA_PAPER/RESCRIPT/SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.qza reference-reads.qza

ln -s /Users/cigom/Documents/MEIOFAUNA_PAPER/RESCRIPT/SSURef_NR99-138.1-tax-derep-uniq.qza reference-tax.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  reference-reads.qza \
  --i-reference-taxonomy reference-tax.qza \
  --o-classifier classifier.qza
```

## 6) Classify
https://btep.ccr.cancer.gov/docs/qiime2/Lesson4/
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads dna-sequences.qza \
  --o-classification dna-sequences-taxonomy.qza  

# Summary res

qiime metadata tabulate \
  --m-input-file dna-sequences-taxonomy.qza \
  --o-visualization dna-sequences-taxonomy.qzv
  
# And import

qiime tools export --input-path dna-sequences-taxonomy.qza --output-path taxonomy_dir

```

## Export qza to text

```bash
# Tax
qiime tools export --input-path SSURef_NR99-138.1-tax-euk.qza --output-path SSURef_NR99-138.1

# 2) Fasta
qiime tools export --input-path SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.qza --output-path SSURef_NR99-138.1
```
Test Domains
```bash
awk '{print $2}' SSURef_NR99-138.1/taxonomy.tsv | sort | uniq -c

20389 d__Archaea;
431329 d__Bacteria;
58790 d__Eukaryota;
```
## Enrich DB
```bash
# 1) SILVA
# 1.1) Tax
qiime tools export --input-path SSURef_NR99-138.1-tax.qza --output-path SSURef_NR99-138.1

# 1.2) Fasta
qiime tools export --input-path SSURef_NR99-138.1-dna-seqs.qza --output-path SSURef_NR99-138.1

# 2) Cat FASTA and TAX
cat SSURef_NR99-138.1.tax pr2_version_4.14.0_SSU_mothur.tax > SSURef_NR99-138.1.pr2_version_4.14.0_SSU.tax

seqkit replace -p "\s.+" SSURef_NR99-138.1-dna-seqs.fasta > SSURef_NR99-138.1-dna-seqs.fa

cat SSURef_NR99-138.1-dna-seqs.fa pr2_version_4.14.0_SSU_mothur.fasta >
SSURef_NR99-138.1.pr2_version_4.14.0_SSU.fasta

# Note:

grep -c "Eukaryota" SSURef_NR99-138.1.pr2_version_4.14.0_SSU.tax # 248,425
wc -l SSURef_NR99-138.1.pr2_version_4.14.0_SSU.tax # 708,111

# 3) Import

fasta_file=SSURef_NR99-138.1.pr2_version_4.14.0_SSU.fasta
tax_file=SSURef_NR99-138.1.pr2_version_4.14.0_SSU.tax

# Sequence
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $fasta_file \
  --output-path ${fasta_file%.fasta}-seqs.qza

# and tax

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $tax_file \
  --output-path ${tax_file%.tax}-tax.qza
  
# Optional 
qiime rescript evaluate-seqs \
    --i-sequences ${fasta_file%.fasta}.qza \
    --o-visualization ${fasta_file%.fasta}.qzv
  
# 4)  
# Redoing from steps 2 and below

sequence=${fasta_file%.fasta}-seqs.qza
taxa=${tax_file%.tax}-tax.qza

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences $sequence \
    --i-taxonomy $taxa \
    --p-labels Eukaryota \
    --p-min-lens 700 \
    --o-filtered-seqs ${sequence%.qza}-euk.qza \
    --o-discarded-seqs ${sequence%.qza}-discarded-seqs-euk.qza
    
# 5) Derep

qiime rescript dereplicate \
    --i-sequences ${sequence%.qza}-euk.qza \
    --i-taxa $taxa \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ${sequence%.qza}-derep.qza \
    --o-dereplicated-taxa ${taxa%.qza}-derep-uniq.qza
    
# 6) Make specific amplicon-region  DB

qiime feature-classifier extract-reads \
    --i-sequences ${sequence%.qza}-derep.qza \
    --p-f-primer TTGTACACACCGCCC \
    --p-r-primer CCTTCYGCAGGTTCACCTAC \
    --p-identity 0.7 \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads ${sequence%.qza}-1389f-1510r.qza
    
# 7) CreateClassify
    
```

# ASV Classification
https://forum.qiime2.org/t/using-q2-clawback-to-assemble-taxonomic-weights/5859/3

```bash
# 4,826 ASVS
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path dna-sequences.fasta \
  --output-path dna-sequences.qza

# 

ln -s /Users/cigom/Documents/MEIOFAUNA_PAPER/RESCRIPT/SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.qza reference-reads.qza

ln -s /Users/cigom/Documents/MEIOFAUNA_PAPER/RESCRIPT/SSURef_NR99-138.1-tax-derep-uniq.qza reference-tax.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  reference-reads.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax-derep-uniq.qza \
  --o-classifier silva-138.1-ssu-nr99-classifier.qza

```