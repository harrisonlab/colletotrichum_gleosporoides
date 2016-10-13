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
  cd $ProjDir
```
