# HiCMCL
High Performance Multi Resolution Markov Clustering for HiC data. Used https://micans.org/mcl/ for the MCL algorithm. 

All the required codes are in the bin folder. The HiCMCL requires python3.4 to run the scripts 

1. BAM to bedpe 

bamtobed tool from bedtools to convert bam to bedpe (bed file with paired end information)

bedtools bamtobed -bedpe -i <BAM>

https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html
  
2. sqlite database dump from bedpe

sql_dump.py

3. Calculating the appropriate parameters (approriate bins and inflation value) for identifying topological boundaries

parametermcl.py

4. Identifying the topological boundaries and classfying them based on their strength

mclclust.py
