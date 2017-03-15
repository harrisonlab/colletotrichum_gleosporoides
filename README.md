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
  # Oxford nanopore 16/12/16
  RawDatDir=/home/miseq_data/minion/2016/minION-Stuart/Colletotrichum/downloaded/pass
  Species=C.gloeosporioides
  Strain=CGMCC3_17371
  mkdir -p raw_dna/minion/$Species/$Strain/16-12-16
  poretools fastq $RawDatDir/ | gzip -cf > raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16.fastq.gz
  poretools stats $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16.stats.txt
  poretools hist $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16.hist
  RawDatDir=/home/miseq_data/minion/2016/minION-Stuart/Colletotrichum/downloaded/fail
  Species=C.gloeosporioides
  Strain=CGMCC3_17371
  mkdir -p raw_dna/minion/$Species/$Strain/16-12-16
  poretools fastq $RawDatDir/ | gzip -cf > raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16_fail.fastq.gz
  poretools stats $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16_fail.stats.txt
  poretools hist $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16_fail.hist
  cat raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16.fastq.gz raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16_fail.fastq.gz > raw_dna/minion/$Species/$Strain/"$Strain"_16-12-16_pass-fail.fastq.gz
  # Oxford nanopore 07/03/17
  RawDatDir=/home/miseq_data/minion/2017/Colletotrichum-2/downloaded/pass
  Species=C.gloeosporioides
  Strain=CGMCC3_17371
  Date=07-03-17
  mkdir -p raw_dna/minion/$Species/$Strain/$Date
  for Fast5Dir in $(ls -d $RawDatDir/*); do
    poretools fastq $Fast5Dir | gzip -cf
  done > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_pass.fastq.gz
  # poretools stats $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.stats.txt
  # poretools hist $RawDatDir/ > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.hist
  # cat raw_dna/minion/$Species/$Strain/"$Strain"_"$Date".fastq.gz raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_fail.fastq.gz > raw_dna/minion/$Species/$Strain/"$Strain"_"$Date"_pass-fail.fastq.gz
```

Pass nanopore data:
This equates to ~ 5.5X coverage, assuming a 57Mb genome
```
  total reads	76476
  total base pairs	312544262
  mean	4086.83
  median	3776
  min	100
  max	23631
  N25	7722
  N50	5861
  N75	4054
```

#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality of illumina reads was visualised using fastqc:
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

<!-- Data quality was also visualised for minion data using fastqc:

```bash
  for RawData in $(ls raw_dna/minion/*/*/*.fastq.gz); do
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData
    qsub $ProgDir/run_fastqc.sh $RawData
  done
``` -->

# Assembly
<!--
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
``` -->


### Canu assembly

```bash
  Organism=C.gloeosporioides
  Strain=CGMCC3_17371
  Run1=$(ls raw_dna/minion/$Organism/$Strain/CGMCC3_17371_16-12-16.fastq.gz)
  Run2=$(ls raw_dna/minion/$Organism/$Strain/CGMCC3_17371_07-03-17_pass.fastq.gz)
  GenomeSz="57m"
  Prefix="$Strain"
  # OutDir="assembly/canu-1.4/$Organism/$Strain"
  OutDir=assembly/canu-1.4/$Organism/"$Strain"_nanopore
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  # qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
  qsub $ProgDir/submit_canu_minion_2lib.sh $Run1 $Run2 $GenomeSz $Prefix $OutDir
```


### Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/canu-1.4/C.gloeosporioides/CGMCC3_17371*/CGMCC3_17371.contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=assembly/canu-1.4/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/canu-1.4/*/*/*.contigs.fasta | grep 'nanopore'); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev | sed 's/_nanopore//g')
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1 | tail -n1);
    TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1 | tail -n1);
    TrimF2_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n2 | tail -n1);
    TrimR2_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n2 | tail -n1);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    InDir=Dir=$(dirname $Assembly)
    OutDir=$InDir/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon_2_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir
  done
```

<!-- ```bash
  Assembly=assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/pilon/pilon.fasta
  # Assembly=assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/scaffolds.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2_spades_pacbio_3
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
``` -->

Inspection of flagged regions didn't identify any contigs that needed to be broken.

### Hybrid assembly:

#### Hybrid assembly: Spades Assembly

```bash
  Organism=C.gloeosporioides
  Strain=CGMCC3_17371
  MinionReads_1=$(ls raw_dna/minion/$Organism/$Strain/*.fastq.gz | grep -v 'fail' | head -n1 | tail -n1)
  MinionReads_2=$(ls raw_dna/minion/$Organism/$Strain/*.fastq.gz | grep -v 'fail' | head -n2 | tail -n1)
  IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
  TrimF1_Read=$(ls $IlluminaDir/F/*.fq.gz);
  TrimR1_Read=$(ls $IlluminaDir/R/*.fq.gz);
  echo $MinionReads_1
  echo $MinionReads_2
  echo $TrimR1_Read
  echo $TrimR1_Read
  OutDir=assembly/spades_minion/$Organism/"$Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
  qsub $ProgDir/subSpades_minion_2lib.sh $MinionReads_1 $MinionReads_2 $TrimF1_Read $TrimR1_Read $OutDir
  # Coverage cuttoff could be set at 73/2 where 73 is the assumed coverage
  # CoverageCuttoff=35
  # OutDir=assembly/spades_pacbio/$Organism/"$Strain"_73x
  # qsub $ProgDir/sub_spades_minion_2lib.sh $MinionReads_1 $MinionReads_2 $TrimF1_Read $TrimR1_Read $OutDir $CoverageCuttoff
```

### Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades_minion/C.gloeosporioides/CGMCC3_17371*/contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=($dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
<!--
Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/canu-1.4/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371.contigs.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/*.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/*.fq.gz);
    OutDir=assembly/spades_pacbio/$Organism/$Strain/pilon
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

```bash
  Assembly=assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/pilon/pilon.fasta
  # Assembly=assembly/spades_pacbio/F.oxysporum_fsp_cepae/Fus2_3/scaffolds.fasta
  Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Fus2_spades_pacbio_3
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```

Inspection of flagged regions didn't identify any contigs that needed to be broken. -->



## Merging pacbio and hybrid assemblies

```bash
  for PacBioAssembly in $(ls assembly/canu-1.4/*/*/polished/*.fasta); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    HybridAssembly=$(ls assembly/spades_*/$Organism/$Strain/contigs.fasta | grep 'minion')
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_spades_first
    AnchorLength=500000
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $HybridAssembly $PacBioAssembly $OutDir $AnchorLength
  done
```

### Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


This merged assembly was polished using Pilon

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta | grep 'spades_first'); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev | cut -f1 -d '_')
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n3 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n3 | tail -n1);
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir/polished
done
```

<!--
Contigs were renamed in accordance with ncbi recomendations

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  # printf "contig_17\tsplit\t780978\t780971\tcanu:missassembly\n"
  for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev | cut -f1 -d '_')
    OutDir=$(dirname $Assembly)
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
``` -->



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

for RawData in $(ls qc_rna/paired/*/*/*/*.fq.gz); do
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

Note this step was run through qlogin

```bash
  qlogin -pe smp 16
```

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'CGMCC3_17371'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
AcceptedHits=$(ls alignment/$Organism/$Strain/accepted_hits.bam)
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
echo "$Organism - $Strain"
mkdir -p $OutDir
cufflinks -o $OutDir/cufflinks -p 16 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
done
```

Output from stdout included:
```
> Processed 50845 loci.                        [*************************] 100%
> Map Properties:
>       Normalized Map Mass: 12686576.36
>       Raw Map Mass: 12686576.36
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 189.73
>                  Estimated Std Dev: 33.02
[15:00:42] Assembling transcripts and estimating abundances.
```

The Estimated Mean: 189.73 allowed calculation of of the mean insert gap to be
5bp 189-(97*2) where 97 was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.


Then Rnaseq data was aligned to each genome assembly:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/*/*); do
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
OutDir=alignment/$Organism/$Strain/$Timepoint
InsertGap='5'
InsertStdDev='33'
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


```bash
# for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep 'Fus2' | grep -e 'Fus2_canu_new' -e 'Fus2_merged' | grep 'cepae' | grep 'Fus2_merged'); do
# for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'proliferatum'); do
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
done
printf "\n"
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
OutDir=gene_pred/braker/$Organism/"$Strain"_braker_new
AcceptedHits=$(ls alignment/$Organism/$Strain/*/accepted_hits.bam)
GeneModelName="$Organism"_"$Strain"_braker_new
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/*/*_braker_new/*/augustus.gff); do
getAnnoFasta.pl $File
OutDir=$(dirname $File)
echo "##gff-version 3" > $OutDir/augustus_extracted.gff
cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls alignment/$Organism/$Strain/*/accepted_hits.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuary:

```bash

for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Jobs=$(qstat | grep 'sub_cuff' | wc -l)
while [ $Jobs -ge 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'sub_cuff' | wc -l)
done
echo "$Organism - $Strain"
OutDir=gene_pred/codingquary/$Organism/$Strain
CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
for BrakerGff in $(ls gene_pred/braker/*/*_braker_new/*/augustus.gff3); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/codingquary/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary

$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

done
```

The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/codingquary/*/*/final); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```



## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP
 * D) From ORF fragments - Signal peptide & RxLR motif  
 * E) From ORF fragments - Hmm evidence of WY domains  
 * F) From ORF fragments - Hmm evidence of RxLR effectors  


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
CurPath=$PWD
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt 20 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```


The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
for SplitDir in $(ls -d gene_pred/final_genes_split/*/*); do
Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
SigpDir=final_genes_spli_signalp-4.1
for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";  
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";  
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";  
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
done
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $TmHeaders
SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
done
```


### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

```bash
for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
cat $File | grep 'Effector' | cut -f1 > $Headers
Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
cat $OutFileHeaders | wc -l
Gff=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended.gff3)
EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
done
```

```
  C.gloeosporioides - CGMCC3_17371
  485
```


#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen output detailing the progress of submission of interproscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
screen -a
cd /home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Genes in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
echo $Genes
$ProgDir/sub_interproscan.sh $Genes
done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteins $InterProRaw
done
```


## B) SwissProt


```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../uniprot/swissprot
SwissDbName=uniprot_sprot
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

```bash
for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl); do
Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
done
```



## 5.2 Identifying PHIbase homologs

The PHIbase database was searched against the assembled genomes using tBLASTx.

```bash
# mkdir -p blast_homology/PHIbase
# cp ../fusarium/analysis/blast_homology/PHIbase/PHI_36_accessions.fa analysis/blast_homology/PHIbase/PHI_36_accessions.fa
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
# Version 4.2 October 3rd 2016
Version=v4.2
DbDir=/home/groups/harrisonlab/phibase/$Version
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/blast_pipe.sh $DbDir/PHI_accessions.fa protein $Assembly
done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness.
<!--
```bash
PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
printf "header\n" > $PhibaseDir/PHI_headers.csv
cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/\r//g' >> $PhibaseDir/PHI_headers.csv
printf "effect\n" > .$PhibaseDir/PHI_virulence.csv
cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | rev | cut -f1 -d '|' | rev  >> $PhibaseDir/PHI_virulence.csv
```


```bash
PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
PhibaseHeaders=$PhibaseDir/PHI_headers.csv
PhibaseVirulence=$PhibaseDir/PHI_virulence.csv
for BlastCSV in $(ls analysis/blast_homology/F*/*/*_PHI_36_accessions.fa_homologs.csv); do
Strain=$(echo $BlastCSV | rev | cut -f2 -d'/' | rev)
echo "$Strain"
OutDir=$(dirname $BlastCSV)
paste -d '\t' $PhibaseHeaders $PhibaseVirulence $BlastCSV | cut -f-3,1185- > $OutDir/"$Strain"_PHIbase_virulence.csv
cat $OutDir/"$Strain"_PHIbase_virulence.csv | grep 'NODE_' | cut -f2 | sort | uniq -c | tee $OutDir/"$Strain"_PHIbase_virulence.txt
done
```
-->


## C) CAZY proteins

Carbohydrte active enzymes were identified using CAZY following recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and  gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
echo "number of CAZY genes identified:"
cat $CazyHeaders | wc -l
# Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
echo "number of Secreted CAZY genes identified:"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done
```

```
number of CAZY genes identified:
1048
530
```

Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)



## Identifying Chittin-masking genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemblies.

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
echo $Assembly
Query=analysis/blast_homology/LysM_genes/Mentlak_et_al_2012_LysM_Sup1.fa
qsub $ProgDir/blast_pipe.sh $Query protein $Assembly
done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
for BlastHits in $(ls analysis/blast_homology/*/*/*Mentlak_et_al_2012_LysM_Sup1.fa_homologs.csv); do
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
Column2=LysM_homolog
NumHits=3
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
```

The M. oryzae Chittin masking gene that has been shown to be essential for
infection was extracted from blast hits:

```bash
Hits=analysis/blast_homology/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371_Mentlak_et_al_2012_LysM_Sup1.fa_homologs.csv
echo "The hit for spl1 is:"
cat $Hits | grep -e 'No.hits' -e 'MGG_10097_Magnaporthe_oryzae'
echo "The hit for spl2 is:"
cat $Hits | grep -e 'No.hits' -e 'MGG_03468_Magnaporthe_oryzae'
```
This identified two hits in the genome, one to NODE 256 and one to node 116.
Node 116 had the better match. The Mentlak paper describes there being two copies
of the Spl gene in the Mo genome, with the second being non-essential for infection.

Hits from spl queries showed greater homology from spl1 (~76% length and 0.49%
identity) than those from spl2 (~47% length and 0.17% identity)

The gff annotations of blast hits were intersected with predicted gene models
to determine which genes represented these regions.

```bash
for HitsGff in $(ls analysis/blast_homology/*/*/*Mentlak_et_al_2012_LysM_Sup1.fa_homologs.gff | grep 'CGMCC3_17371'); do
Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended.gff3)
OutDir=$(dirname $HitsGff)
LysmIntersect=$OutDir/"$Strain"_LysM_hit_genes.bed
bedtools intersect -wao -a $HitsGff -b $GeneGff > $LysmIntersect
echo "Gene models intersected are:"
cat $LysmIntersect | grep -w 'gene' | cut -f18 | cut -f2 -d '=' | tr -d ';' | sort | uniq
done
```

Gene models intersected were g11852 & g14420.
```
>g11852.t1
MQTSYIFTTLLAAAGLVAALPQATPTQVTPTGTASATASATPTCSQGPVVDYTVVSGDTL
TIISQKLSSGICDIAKLNSLENPNLILLGQVLKVPTHICNPDNTSCLSKPGTATCVEGGP
ATYTIQKGDTFFIVAGDLGLDVNALLAANEGVDPLLLQEGQVINIPVCKX
>g14420.t1
MQFSIFTVLAAAASAVVALPAATPTAAATATPSATCGKIGNFHKTTVKAGQTLTTIAQRY
NSGICDIAWQNKLANPNVIFAGQVLLVPVDVCNPDNTSCITPTGEATCVTGGPATYTIKS
GDTFFVVAQSLGITTDSLTGANPGVAAESLQVGQVIKVPVCX
```

These results were cross referenced against interproscan annotation results:

```bash
  OutDir=analysis/LysM
  mkdir -p $OutDir
  InterproTsv=gene_pred/interproscan/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371_interproscan.tsv
  cat $InterproTsv | grep -w 'LysM domain' | cut -f1 | sort | uniq > $OutDir/LysM_gene_headers.txt
  echo "Number of LysM proteins:"
  cat $OutDir/LysM_gene_headers.txt | wc -l
  Secretedheaders=gene_pred/final_genes_signalp-4.1/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371_final_sp_no_trans_mem_headers.txt
  cat $Secretedheaders | grep -f $OutDir/LysM_gene_headers.txt > $OutDir/LysM_gene_headers_secreted.txt
  echo "Number of secreted LysM proteins:"
  cat $OutDir/LysM_gene_headers_secreted.txt | wc -l
```

Number of LysM proteins:
31
Number of secreted LysM proteins:
18

The secreted list included both g11852 & g14420.

## Extracting sequence data for phylogenetic loci

Nucleotide sequences of phylogenetic loci ACT, CAL, GAPD, ITS, GS, Btub, ApMat
were identified in the genome assemblies.

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
echo $Assembly
Query=analysis/blast_homology/phylogenetic_loci/phylogenetic_loci.fa
qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
for BlastHits in $(ls analysis/blast_homology/*/*/*phylogenetic_loci.fa_homologs.csv); do
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
Column2=phylogenetic_locus
NumHits=4
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
```

The hits of these loci were manually edited so that the gff annotation name
contained the locus ID

```bash
  # for CGMCC3_17371
  HitsGff=analysis/blast_homology/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371_phylogenetic_loci.fa_homologs.gff
  EditedHitsGff=analysis/blast_homology/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371_phylogenetic_loci_edited_homologs.gff
  cp  $HitsGff $EditedHitsGff
  nano $EditedHitsGff
  # for nara_gc5
  HitsGff=analysis/blast_homology/C.gloeosporioides/Nara_gc5/Nara_gc5_phylogenetic_loci.fa_homologs.gff
  EditedHitsGff=analysis/blast_homology/C.gloeosporioides/Nara_gc5/Nara_gc5_phylogenetic_loci_edited_homologs.gff
  cp  $HitsGff $EditedHitsGff
  nano $EditedHitsGff
```

The hit and the 500bp upstream and downstream
 of the hit were extracted as a fasta file:

```bash
  for HitsGff in $(ls analysis/blast_homology/C.gloeosporioides/*/*_phylogenetic_loci_edited_homologs.gff); do
    Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Genome=$(ls repeat_masked/$Organism/$Strain/first_assembly/"$Strain"_contigs_unmasked.fa)
    # cat $Genome | sed 's/gi|.*|gb|A//g' | sed 's/ Colletotrichum gloeosporioides.*//g' > $OutDir/genome_parsed.fa
    OutDir=analysis/phylogenetics/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
    $ProgDir/gffexpander.pl +- 500 $HitsGff > $OutDir/"$Strain"_loci_exp_500bp.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/sequence_extractor.py --coordinates_file $OutDir/"$Strain"_loci_exp_500bp.gff --header_column 1 --start_column 4 --stop_column 5 --strand_column  7 --id_column 9 --fasta_file $Genome | sed 's/"ID=//g' | sed 's/";//g'> $OutDir/"$Strain"_loci_exp_500bp.fa
    # $ProgDir/sequence_extractor.py --coordinates_file $OutDir/"$Strain"_loci_exp_500bp.gff --header_column 1 --start_column 4 --stop_column 5 --strand_column  7 --id_column 9 --fasta_file $OutDir/genome_parsed.fa | sed 's/"ID=//g' | sed 's/";//g'> $OutDir/"$Strain"_loci_exp_500bp.fa
  done
```
