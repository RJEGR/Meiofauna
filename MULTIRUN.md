Ricardo Gomez Reyes
Preprocessing step
```bash
for f in $(ls *R1* | grep fastq); do ./trim_galore -q 10 --paired ${f%_R*}_R1_001.fastq.gz ${f%_R*}_R2_001.fastq.gz; done
```

**1) Importing data to qiime**
See qiime fastq manifest [first](https://docs.qiime2.org/2024.2/tutorials/importing/#fastq-manifest-formats) prior to import data. 
```bash
for i in $(ls *R1*fq.gz)
do
withpath="${i}"
fname=${withpath##*/}
bs="${fname%_R*}"
echo "${bs%_*}\t" `printf "$PWD/${bs}_R1_001_val_1.fq.gz\t$PWD/${bs}_R2_001_val_2.fq.gz"`
done > manifest.tsv

# vi manifest and paste header line
# sample-id     forward-absolute-filepath       reverse-absolute-filepath

```
Using the Phred 33 offset as Hiseq-2500 series use Illumina 1.8 or later software. Review more quality score variant [here](https://scikit.bio/docs/dev/generated/skbio.io.format.fastq.html#quality-score-variants) prior to run the follow code. 
```bash

# scp -r acisterna@ron.sr.unh.edu://home/gomre/acisterna/qiime_analyses_Alejandro/qiime2_alejandro/analisis-informe-final-cigom/raw-seqs/TRIMMED/LIBRARY

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path MANIFEST.tsv \
  --output-path demultiplexed-sequences.qza
```

summarize
```bash
qiime demux summarize \
  --i-data demultiplexed-sequences.qza \
  --o-visualization demultiplexed-sequences-summ.qzv
```

**2) Denoising sequence data with DADA2**
```bash
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demuxsequences.qza \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --p-n-threads 2 \
  --o-representative-sequences asv-sequences.qza \
  --o-table feature-table.qza \
  --o-denoising-stats dada2-stats.qza  

```