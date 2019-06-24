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

<!--
The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/final/*/*/final | grep -v '_braker'); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```

```
gene_pred/final/C.gloeosporioides/CGMCC3_17371/final
16758
1689
18447
``` -->

<!--
In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.


The next step had problems with the masked pacbio genome. Bioperl could not read in
the fasta sequences. This was overcome by wrapping the unmasked genome and using this
fasta file.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
NewName=$(echo $Assembly | sed 's/_unmasked.fa/_unmasked_wrapped.fa/g')
cat $Assembly | fold > $NewName
done
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep -v '_braker'); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev | sed 's/_UTR//g')
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
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_unmasked_wrapped.fa)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta
done
```

Duplicated transcript:CUFF_3819_2_309.t1

```bash
  for Gff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3); do
  	Strain=$(echo $Gff | rev | cut -d '/' -f3 | rev)
  	Organism=$(echo $Gff | rev | cut -d '/' -f4 | rev)
  	echo "$Strain - $Organism"
  	cat $Gff | grep -w 'gene' | wc -l
  done
```

```
CGMCC3_17371 - C.gloeosporioides
18148
```

The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/final/*/*/final | grep -v '_braker'); do
echo $DirPath;
cat $DirPath/final_genes_appended_renamed.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_appended_renamed.pep.fasta | grep '>' | grep '.t1' | wc -l;
echo "";
done
```

```
gene_pred/final/C.gloeosporioides/CGMCC3_17371/final
18446
18147
``` -->

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
