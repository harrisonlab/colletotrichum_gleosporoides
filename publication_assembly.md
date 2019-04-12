# Colletotrichum gloeosporioides
==========

Scripts used for the analysis of Colletotrichum gloeosporioides genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides


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
```

## Estimating sequencing coverage

For Minion data:
```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=57
OutDir=$(dirname $RawData)
mkdir -p $OutDir
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*_cov.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

MinION coverage

```
CGMCC3_17371	13.76
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



Read coverage was estimated from the trimemd datasets:

```bash
GenomeSz=55
for Reads in $(ls qc_dna/paired/C.gloeosporioides/CGMCC3_17371/*/*.fq.gz); do
echo $Reads
OutDir=$(dirname $Reads)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $Reads $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/paired/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*/*_cov.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```bash
  CGMCC3_17371	72.98
```


Splitting reads and trimming adapters using porechop
```bash
for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz | grep -e '07-03-17_pass.fastq' -e '16-12-16.fastq'); do
Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RawReads | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=qc_dna/minion/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
done
```


The two datasets were combined into a single dataset for later assembly correction:

```bash
  Reads1=$(ls qc_dna/minion/C.gloeosporioides/CGMCC3_17371/*.fastq.gz | grep -v 'appended' | head -n1 | tail -n1)
  Reads2=$(ls qc_dna/minion/C.gloeosporioides/CGMCC3_17371/*.fastq.gz | grep -v 'appended' | head -n2 | tail -n1)
  OutDir=$(dirname $Reads1)
  cat $Reads1 $Reads2 > $OutDir/minion_appended.fastq.gz
```


# Assembly


### Canu assembly

```bash
  Organism=C.gloeosporioides
  Strain=CGMCC3_17371
  Reads1=$(ls qc_dna/minion/$Organism/$Strain/*.fastq.gz | grep -v 'appended' | head -n1 | tail -n1)
  Reads2=$(ls qc_dna/minion/$Organism/$Strain/*.fastq.gz | grep -v 'appended' | head -n2 | tail -n1)
  GenomeSz="57m"
  Prefix="$Strain"
  OutDir=assembly/canu-1.8/$Organism/"$Strain"
  mkdir -p $OutDir
  AdditionalCommands=$PWD/$OutDir/additional_commands.txt
  cat "correctedErrorRate 0.16" > $AdditionalCommands
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu_minion_2lib.sh $Reads1 $Reads2 $GenomeSz $Prefix $OutDir $AdditionalCommands
```

Quast and busco were run to assess assembly quality:

```bash
for Assembly in $(ls assembly/canu-1.8/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371.contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
# Quast
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.8/*/*/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```

### Hybrid assembly:

#### Hybrid assembly: Spades Assembly

```bash
  Organism=C.gloeosporioides
  Strain=CGMCC3_17371
  MinionReads_1=$(ls qc_dna/minion/$Organism/$Strain/*.fastq.gz | grep -v 'fail' | head -n1 | tail -n1)
  MinionReads_2=$(ls qc_dna/minion/$Organism/$Strain/*.fastq.gz | grep -v 'fail' | head -n2 | tail -n1)
  IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
  TrimF1_Read=$(ls $IlluminaDir/F/*.fq.gz);
  TrimR1_Read=$(ls $IlluminaDir/R/*.fq.gz);
  echo $MinionReads_1
  echo $MinionReads_2
  echo $TrimR1_Read
  echo $TrimR1_Read
  OutDir=assembly/spades-3.11_minion/$Organism/"$Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
  qsub $ProgDir/subSpades_minion_2lib.sh $MinionReads_1 $MinionReads_2 $TrimF1_Read $TrimR1_Read $OutDir
```


Quast and busco were run to assess assembly quality:

```bash
for Assembly in $(ls assembly/spades-3.11_minion/C.gloeosporioides/CGMCC3_17371*/contigs.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
# Quast
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
