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

Note - Genome size was determined by Assembly size of previous Cg sequencing project:
https://www.ncbi.nlm.nih.gov/genome?LinkName=nuccore_genome&from_uid=530480622


```bash
Size=55
for TrimPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
TrimF=$(ls $TrimPath/F/*.fq.gz)
TrimR=$(ls $TrimPath/R/*.fq.gz)
echo $TrimF
echo $TrimR
mkdir -p tmp_coverage
cat $TrimF | gunzip -cf > tmp_coverage/F_read.fq
cat $TrimR | gunzip -cf > tmp_coverage/R_read.fq
count_nucl.pl -i tmp_coverage/F_read.fq -i tmp_coverage/R_read.fq -g $Size
rm -r tmp_coverage
done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

This was performed for strains with single runs of data

```bash
for TrimPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
TrimF=$(ls $TrimPath/F/*.fq.gz)
TrimR=$(ls $TrimPath/R/*.fq.gz)
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

```bash
for StrainPath in $(ls -d qc_dna/paired/*/*); do
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
	Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
	F_Read=$(ls $StrainPath/F/*.fq.gz)
	R_Read=$(ls $StrainPath/R/*.fq.gz)
	OutDir=assembly/spades/$Organism/$Strain
	echo $F_Read
	echo $R_Read
	qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 15
done
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $BestAss | rev | cut -f3 -d '/' | rev)
Organism=$(echo $BestAss | rev | cut -f4 -d '/' | rev)
OutDir=repeat_masked/$Organism/"$Strain"/first_assembly
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and harmasked files.

```bash
for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/*/*/first_assembly); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
```
```
C.gloeosporioides	CGMCC3_17371
The number of bases masked by RepeatMasker:	1226490
The number of bases masked by TransposonPSI:	269923
The total number of masked bases are:	1411323

C.gloeosporioides	Nara_gc5
The number of bases masked by RepeatMasker:	756936
The number of bases masked by TransposonPSI:	195889
The total number of masked bases are:	933627
```


# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

# Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
for Genome in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
echo $Genome;
qsub $ProgDir/sub_cegma.sh $Genome dna;
done
```

Outputs were summarised using the commands:
```bash
for File in $(ls gene_pred/cegma/*/*/*_dna_cegma.completeness_report); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
Species=$(echo $File | rev | cut -f3 -d '/' | rev);
printf "$Species\t$Strain\n";
cat $File | head -n18 | tail -n+4;printf "\n";
done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```


#Gene prediction

Gene prediction was performed for Fusarium genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* qc of RNA seq data is detailed below:

```bash
  # RNAseq data from Zhang et al 2015 was used to annotate the genome
  Species=C.gloeosporioides
  Strain=ES026
  OutDir=raw_rna/paired/$Species/$Strain
  mkdir -p $OutDir/paired
  fastq-dump -O $OutDir/paired --split-files SRR1171641
  # fastq-dump -O $OutDir/paired --gzip --split-files SRR1171641
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cat $OutDir/paired/SRR1171641_1.fastq | gzip -cf > $OutDir/F/SRR1171641_1.fastq.gz
  cat $OutDir/paired/SRR1171641_2.fastq | gzip -cf > $OutDir/R/SRR1171641_2.fastq.gz
```


Perform qc of RNAseq timecourse data
```bash
for FilePath in $(ls -d raw_rna/paired/*/*); do
echo $FilePath
FileF=$(ls $FilePath/F/*.gz)
FileR=$(ls $FilePath/R/*.gz)
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
done
```

Data quality was visualised using fastqc:
```bash

	for RawData in $(ls qc_rna/paired/F.oxysporum_fsp_cepae/*/*/*.fq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

#### Aligning

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

```bash
	# for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
# Jobs=$(qstat | grep 'rna_qc_fas' | wc -l)
# while [ $Jobs -ge 1 ]; do
#   sleep 5m
#   printf "."
#   Jobs=$(qstat | grep 'rna_qc_fas' | wc -l)
# done
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/*/*); do
  # Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
  # echo "$Timepoint"
  FileF=$(ls $RNADir/F/*_trim.fq.gz)
  FileR=$(ls $RNADir/R/*_trim.fq.gz)
  # OutDir=alignment/$Organism/$Strain/$Timepoint
  OutDir=alignment/$Organism/$Strain
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
  qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
done
done
```
Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
AcceptedHits=$(ls alignment/$Organism/$Strain/*/concatenated.bam)
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
echo "$Organism - $Strain"
mkdir -p $OutDir
cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
done
```

Output from stdout included:
```
	Processed 22484 loci.                        [*************************] 100%
	Map Properties:
	     Normalized Map Mass: 50507412.55
	     Raw Map Mass: 50507412.55
	     Fragment Length Distribution: Empirical (learned)
	                   Estimated Mean: 181.98
	                Estimated Std Dev: 78.39
	[13:02:48] Assembling transcripts and estimating abundances.
	Processed 22506 loci.                        [*************************] 100%
```

The Estimated Mean: 181.98 allowed calculation of of the mean insert gap to be
-20bp 182-(97*2) where 97 was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.

<!--
Then Rnaseq data was aligned to each genome assembly:

```bash
# for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'ncbi' | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'A8'); do
# for Assembly in $(ls assembly/merged_canu_spades/*/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta); do
for Assembly in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -e 'FO47_72hrs_rep2'); do
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
OutDir=alignment/$Organism/$Strain/$Timepoint
InsertGap='-20'
InsertStdDev='78'
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
done
printf "\n"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
done
done
```
 -->
