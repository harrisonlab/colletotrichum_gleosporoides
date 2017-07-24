# Reference genomes

This file details the commands to download and analyse reference genomes for
Colletotrichum gloeosporioides and related species.



#Building of directory structure

```bash
  ProjDir=/home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides
  cd $ProjDir

  # Data from https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=ANPB01&display=contigs&page=1
  Species=C.gloeosporioides
  Strain=Nara_gc5
  Source=NCBI
  mkdir -p assembly/external_group/$Species/$Strain/$Source
  cd assembly/external_group/$Species/$Strain/$Source
  wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AN/PB/ANPB01/ANPB01.1.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AN/PB/ANPB01/ANPB01.1.fsa_nt.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AN/PB/ANPB01/ANPB01.1.bbs.gz
  gunzip *.gz
```
<!--
```bash
  ProjDir=/home/groups/harrisonlab/project_files/colletotrichum_gloeosporioides
  cd $ProjDir
  # Data from https://www.ncbi.nlm.nih.gov/assembly/GCA_000149035.1
  Species=C.gramminicola
  Strain=M1_001
  Source=Broad
  mkdir -p assembly/external_group/$Species/$Strain/$Source
  cd assembly/external_group/$Species/$Strain/$Source
  # wget -r ftp://archive.broadinstitute.org/ftp/pub/annotation/fungi/colletotrichum/genomes/colletotrichum_graminicola_m1/
  wget -r http://archive.broadinstitute.org/ftp/pub/annotation/fungi/colletotrichum/genomes/colletotrichum_graminicola_m1/*
  gunzip *.gz
  cd $ProjDir
  Species=C.higgsisnum
  Strain=imi_349063
  Source=Broad
  mkdir -p assembly/external_group/$Species/$Strain/$Source
  cd assembly/external_group/$Species/$Strain/$Source
  wget -r http://archive.broadinstitute.org/ftp/pub/annotation/fungi/colletotrichum/genomes/colletotrichum_higginsianum_imi_349063
  gunzip *.gz
  cd $ProjDir
``` -->

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  for BestAss in $(ls assembly/external_group/*/*/*/ANPB01.1.fsa_nt); do
    Strain=$(echo $BestAss | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $BestAss | rev | cut -f4 -d '/' | rev)
    OutDir=repeat_masked/$Organism/"$Strain"/first_assembly
    qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
    qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
  done
```


```bash
for Assembly in  $(ls assembly/external_group/*/*/*/ANPB01.1.fsa_nt); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/assembly
  qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```
