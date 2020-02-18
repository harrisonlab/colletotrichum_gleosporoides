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

<!--
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
-->

While logged in to the new cluster:

```bash
WorkDir=$HOME/tmp/smartdenovo_Cg


# ---------------
# Step 1
# Collect inputs
# ---------------

ProjDir=/oldhpc/home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides
FastaIn=$(ls $ProjDir/assembly/canu-1.8/*/*/*.trimmedReads.fasta.gz)
Organism=$(echo $FastaIn | rev | cut -f3 -d '/' | rev)
Strain=$(echo $FastaIn | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=$ProjDir/assembly/SMARTdenovo/$Organism/"$Strain"
echo  "Running SMARTdenovo with the following inputs:"
echo "FastaIn - $FastaIn"
echo "Prefix - $Prefix"
echo "OutDir - $OutDir"

# ---------------
# Step 2
# Run SMARTdenovo
# ---------------

mkdir -p $WorkDir
cd $WorkDir
Fasta=$(basename $FastaIn)
cp $FastaIn $Fasta

cat $Fasta | gunzip -cf > reads.fa
smartdenovo.pl -t 40 reads.fa -p $Prefix > $Prefix.mak

make -f $Prefix.mak 2>&1 | tee "$Prefix"_run_log.txt

rm $Prefix.mak
rm $Fasta
rm reads.fa
mkdir -p $OutDir
for File in $(ls $WorkDir/wtasm*); do
  NewName=$(echo $File | sed "s/wtasm/$Prefix/g")
  mv $File $NewName
done
cp -r $WorkDir/* $OutDir/.
rm -r $WorkDir
```

Quast and busco were run to assess assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/C.gloeosporioides/CGMCC3_17371/CGMCC3_17371_smartdenovo.dmo.lay.utg); do
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

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt | grep -e 'canu-1.8' -e '_contigs.txt'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
  short_summary_CGMCC3_17371.contigs.txt	325	0	337	3063	3725
  short_summary_contigs.txt	3685	6	21	19	3725
```

## KAT kmer spectra analysis

Kat was run on the hybrid assembly:

```bash
  for Assembly in $(ls assembly/spades-3.11_minion/C.gloeosporioides/CGMCC3_17371*/contigs.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    ReadsF=$(ls qc_dna/paired/$Organism/$Strain/F/*_trim.fq.gz)
    ReadsR=$(ls qc_dna/paired/$Organism/$Strain/R/*_trim.fq.gz)
    OutDir=$(dirname $Assembly)/kat
    Prefix="repeat_masked"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/kat
    qsub $ProgDir/sub_kat.sh $Assembly $ReadsF $ReadsR $OutDir $Prefix 250
  done
```


### Pilon assembly correction

Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/spades-3.11_minion/C.gloeosporioides/CGMCC3_17371*/contigs.fasta); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/*.fq.gz);
OutDir=$(dirname $Assembly)/pilon
Iterations=5
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

Contigs were renamed
```bash
echo "" > tmp.txt
for Assembly in $(ls assembly/spades-3.11_minion/*/*/pilon/*.fasta | grep 'pilon_5'); do
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess assembly quality:

```bash
for Assembly in $(ls assembly/spades-3.11_minion/C.gloeosporioides/CGMCC3_17371*/pilon/pilon_min_500bp_renamed.fasta); do
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

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt | grep -e 'pilon_min_500bp_renamed'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
short_summary_pilon_min_500bp_renamed.txt	3686	6	21	18	3725
```

## KAT kmer spectra analysis

Kat was run on the hybrid assembly:

```bash
  for Assembly in $(ls assembly/spades-3.11_minion/C.gloeosporioides/CGMCC3_17371*/pilon/pilon_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    ReadsF=$(ls qc_dna/paired/$Organism/$Strain/F/*_trim.fq.gz)
    ReadsR=$(ls qc_dna/paired/$Organism/$Strain/R/*_trim.fq.gz)
    OutDir=$(dirname $Assembly)/kat
    Prefix="repeat_masked"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/kat
    qsub $ProgDir/sub_kat.sh $Assembly $ReadsF $ReadsR $OutDir $Prefix 250
  done
```



A Bioproject and Biosample was made with NCBI genbank for submission of genomes.
Following the creation of these submissions, the .fasta assembly was uploaded
through the submission portal. A note was provided requesting that the assembly
be run through the contamination screen to aid a more detailed resubmission in
future. The returned FCSreport.txt was downloaded from the NCBI webportal and
used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following loactions:

```bash
  for Assembly in $(ls assembly/spades-3.11_minion/C.gloeosporioides/CGMCC3_17371*/pilon/pilon_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
  done
```


These downloaded files were used to correct assemblies:

```bash
for Assembly in $(ls assembly/spades-3.11_minion/C.gloeosporioides/CGMCC3_17371*/pilon/pilon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/Contamination*.txt)
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/report.txt); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    echo;
    echo $Organism;
    echo $Strain;
    cat $Assembly;
  done > assembly/quast_results.txt
```

```
Assembly                   contigs_min_500bp_renamed
# contigs (>= 0 bp)        447                      
# contigs (>= 1000 bp)     308                      
Total length (>= 0 bp)     58056435                 
Total length (>= 1000 bp)  57960885                 
# contigs                  447                      
Largest contig             2278540                  
Total length               58056435                 
GC (%)                     53.20                    
N50                        905947                   
N75                        390361                   
L50                        23                       
L75                        47                       
# N's per 100 kbp          0.00
```


# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

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

```
repeat_masked/C.gloeosporioides/CGMCC3_17371/filtered_contigs/CGMCC3_17371_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
2109446
```


# Gene Prediction


Gene prediction followed three steps:
Pre-gene prediction
- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene model training
- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
Gene prediction
- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

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
  Strain=CGMCC3_17371
  OutDir=/data/scratch/armita/colletotrichum_gloeosporioides/raw_rna/paired/$Species/$Strain
  mkdir -p $OutDir/paired
  fastq-dump -O $OutDir/paired --gzip --split-files SRR5194995
  fastq-dump -O $OutDir/paired --gzip --split-files SRR5194994
  fastq-dump -O $OutDir/paired --gzip --split-files SRR5194993
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  mv $OutDir/paired/*_1.fastq.gz $OutDir/F/.
  mv $OutDir/paired/*_2.fastq.gz $OutDir/R/.
```


Perform qc of RNAseq timecourse data
```bash
cd /data/scratch/armita/colletotrichum_gloeosporioides
for FileF in $(ls $FilePath/F/*.gz); do
echo $FileF
FileR=$(echo $FileF | sed "s/F/R/g" | sed "s/_1.fastq/_2.fastq/g")
ls $FileR
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
done
```

Data quality was visualised using fastqc:
```bash
cd /data/scratch/armita/colletotrichum_gloeosporioides
for RawData in $(ls qc_rna/paired/*/*/*/*.fq.gz); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
echo $RawData;
qsub $ProgDir/run_fastqc.sh $RawData
done
```

#### Aligning

RNAseq data was aligned to the assemblies using STAR to provide evidence for gene models

```bash
  cd /data/scratch/armita/colletotrichum_gloeosporioides

  for Assembly in $(ls ../../../../home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides/repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
  # for Assembly in $(ls ../../../../home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides/assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d qc_rna/paired/*/*); do
      FileNum=$(ls $RNADir/F/*_trim.fq.gz | wc -l)
      for num in $(seq 1 $FileNum); do
        while [ $Jobs -gt 1 ]; do
          sleep 1m
          printf "."
          Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        FileF=$(ls $RNADir/F/*_trim.fq.gz | head -n $num | tail -n1)
        FileR=$(ls $RNADir/R/*_trim.fq.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed "s/_R.*_trim.fq.gz//g")
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        # Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
        Timepoint=$num
        echo "$Timepoint"
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
      done
    done
  done
```


Alignments were concatenated prior to gene prediction. This step was done through
a qlogin session on the cluster.

```bash
	for StrainDir in $(ls -d alignment/star/*/*); do
	BamFiles=$(ls $StrainDir/*/*/star_aligmentAligned.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
	OutDir=$StrainDir/concatenated
	mkdir -p $OutDir
	samtools merge -f $OutDir/all_reads_concatenated.bam $BamFiles
	done
```


#### Braker prediction

Before braker predictiction was performed, I double checked that I had the
genemark key in my user area and copied it over from the genemark install
directory:

```bash
	ls ~/.gm_key
  cp /home/armita/prog/genemark/2019/gm_key_64 ~/.gm_key
```

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
AcceptedHits=$(ls ../../../../../data/scratch/armita/colletotrichum_gloeosporioides/alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
GeneModelName="$Organism"_"$Strain"_braker_new
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3); do
getAnnoFasta.pl $File
OutDir=$(dirname $File)
echo "##gff-version 3" > $OutDir/augustus_extracted.gff
cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls ../../../../../data/scratch/armita/colletotrichum_gloeosporioides/alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```


Secondly, genes were predicted using CodingQuary:

```bash

for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
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

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

Note - Ensure that the "TPSI_appended.fa" assembly file is correct.

```bash
for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
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
# -
# This section is edited
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
$ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
# -
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
done
```

```bash
  for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep -v 'braker'); do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    Gene=$(cat $GffAppended | grep -w 'gene' | wc -l)
    Protein=$(cat $GffAppended | grep -w 'mRNA' | wc -l)
    Augustus=$(cat $GffAppended | grep -w 'gene' | grep 'AUGUSTUS' | wc -l)
    CodingQuary=$(cat $GffAppended | grep -w 'gene' | grep 'CodingQuarry_v2.0' | wc -l)
    printf "$Organism\t$Strain\t$Gene\t$Protein\t$Augustus\t$CodingQuary\n"
  done
```

```
C.gloeosporioides       CGMCC3_17371    18143   18447   16459   1684
```


In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified

<!-- Codingquary was noted to predict a gene that went beyond the end of contig 47 in
isolate 1177.

As such this gene was removed manually:

```bash
GffAppended=$(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '1177')
cp $GffAppended tmp.gff
cat tmp.gff | grep -v 'CUFF_8208_1_74' > $GffAppended
``` -->


```bash
  for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep -v 'braker'); do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/final/$Organism/$Strain/final
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should
    # be changed
    sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
  done
```

Gene CUFF_3819_1_308.t2 was identified as duplicated

```bash
  for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3 | grep -v 'braker'); do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    Gene=$(cat $GffAppended | grep -w 'gene' | wc -l)
    Protein=$(cat $GffAppended | grep -w 'mRNA' | wc -l)
    Augustus=$(cat $GffAppended | grep -w 'gene' | grep 'AUGUSTUS' | wc -l)
    CodingQuary=$(cat $GffAppended | grep -w 'gene' | grep 'CodingQuarry_v2.0' | wc -l)
    printf "$Organism\t$Strain\t$Gene\t$Protein\t$Augustus\t$CodingQuary\n"
  done
```

```
C.gloeosporioides       CGMCC3_17371    18143   18446   16459   1684
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
for Genes in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -v '_braker'); do
echo $Genes
$ProgDir/sub_interproscan.sh $Genes
done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -v '_braker'); do
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
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -v '_braker'); do
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
for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_vMar2018_10_hits.tbl); do
Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_vMar2018_tophit_parsed.tbl
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
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
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -v '_braker'); do
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
SigpDir=final_genes_signalp-4.1
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
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -v '_braker'); do
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

```
C.gloeosporioides - CGMCC3_17371
2070
```


### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -v '_braker'); do
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
612
```


## C) CAZY proteins

Carbohydrte active enzymes were identified using CAZY following recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -v '_braker'); do
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
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
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
