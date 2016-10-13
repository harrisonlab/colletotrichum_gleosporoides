# Colletotrichum gloeosporioides
==========

Scripts used for the analysis of Colletotrichum gloeosporioides genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

<!--
Analyses performed on these genomes involved BLAST searching for:

 ls contigs were identified using:
Alignment of raw reads to assembled genomes
Assembly of remaining reads
-->


#Building of directory structure

```bash
  ProjDir=/home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides
  mkdir -p $ProjDir
  cd $ProjDir
  # Data run on 13/10/16 161010_M04465_0026_000000000-APP5B
  RawDatDir=/home/miseq_readonly/miseq_data/2016/RAW/161010_M04465_0026_000000000-APP5B/Data/Intensities/BaseCalls
  Species=C.gloeosporioides
  Strain=CGMCC3_17371
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R

  cp $RawDatDir/Colletotrichum_S1_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/.
  cp $RawDatDir/Colletotrichum_S1_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/.
```


#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```

Trimming was first performed on all strains that had a single run of data:

```bash
  for StrainPath in $(ls -d raw_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ReadsF=$(ls $StrainPath/F/*.fastq*)
    ReadsR=$(ls $StrainPath/R/*.fastq*)
    echo $ReadsF
    echo $ReadsR
    qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
  done
```

Data quality was visualised once again following trimming:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

The sequencing coverage for isolates was estimated by counting all the
nucleotides in the trimmed fastq files and dividing this by the estimated
genome size.

```bash
  Size=50
  for TrimPath in $(ls -d qc_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    TrimF=$(ls $TrimPath/F/*.fastq*)
    TrimR=$(ls $TrimPath/R/*.fastq*)
    echo $TrimF
    echo $TrimR
    mkdir -p tmp_coverage
    cat $TrimF | gunzip -cf > tmp_coverage/F_read.fq
    cat $TrimR | gunzip -cf > tmp_coverage/R_read.fq
    count_nucl.pl -i tmp_coverage/F_read.fq tmp_coverage/R_read.fq -g $Size
    rm -r tmp_coverage
  done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

This was performed for strains with single runs of data

```bash
	for TrimPath in $(ls -d qc_dna/paired/*/*); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		echo $TrimF
		echo $TrimR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
	done
```


mode kmer abundance prior to error correction was reported using the following
commands:

```bash
for File in $(ls qc_dna/kmc/*/*/*_true_kmer_summary.txt); do
basename $File;
tail -n3 $File | head -n1 ;
done
```
