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

Starting to ...
https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494

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