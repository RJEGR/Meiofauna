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

### SILVA
The SILVA 138 SSU Ref NR 99 database has been downloaded using rescript plugging for subsequent data analysis (*SSURef_Nr99_tax_silva.fasta),

```bash
WD=/Users/cigom/Documents/MEIOFAUNA_PAPER/RESCRIPT

cd $WD

conda activate qiime2-2023.2
qiime rescript --help

# 1) Download db

target="SSURef_NR99"
version="138.1"

qiime rescript get-silva-data \
    --p-version $version \
    --p-target $target \
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

# wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*" # 

# 18S ribosomal RNA sequences (SSU) from Fungi type and reference material (BioProject PRJNA39195) [blast db files <makeblast>]
18S_fungal_sequences.tar.gz
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

## FeatureData conversion (rna to dna)
If you'd like to be able to jump to steps that only take FeatureData[Sequence] as input you can convert your data to FeatureData[Sequence] like so:

```bash
qiime rescript reverse-transcribe \
    --i-rna-sequences SSURef_NR99-138.1-rna-seqs.qza \
    --o-dna-sequences SSURef_NR99-138.1-dna-seqs.qza
```

## Parse taxonomy
```bash
qiime rescript parse-silva-taxonomy \
    --i-taxonomy-ranks SSURef_NR99-138.1-tax.qza \
    --p-no-rank-propagation \
    --p-ranks domain kingdom phylum class order suborder family genus \
    --p-include-species-labels \
    --o-taxonomy SSURef_NR99-138.1-tax-ranks.qza
```

Filter DB by Domain
```bash
# Filtering sequences by length and taxonomy
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences SSURef_NR99-138.1-rna-seqs.qza \
    --i-taxonomy SSURef_NR99-138.1-tax.qza \
    --p-labels Eukaryota \
    --o-filtered-seqs SSURef_NR99-138.1-rna-euk-seqs.qza \
    --o-discarded-seqs SSURef_NR99-138.1-rna-euk-tax.qza 
```

## Dereplicate sequences
```bash
qiime rescript dereplicate \
    --i-sequences silva-138-ssu-nr99-seqs-filt.qza \
    --i-taxa silva-138-ssu-nr99-tax.qza \
    --p-perc-identity 0.97 \
```

## Export qza to text

```bash
# Tax
qiime tools export --input-path SSURef_NR99-138.1-tax.qza --output-path SSURef_NR99-138.1

# 2) Fasta
qiime tools export --input-path SSURef_NR99-138.1-dna-seqs.qza --output-path SSURef_NR99-138.1
```
Test Domains
```bash
awk '{print $2}' SSURef_NR99-138.1-tax/taxonomy.tsv | sort | uniq -c

20389 d__Archaea;
431329 d__Bacteria;
58790 d__Eukaryota;
```