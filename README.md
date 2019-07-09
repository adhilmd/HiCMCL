# HiCMCL
High Performance Multi Resolution Markov Clustering for HiC data

All the required codes are in the bin folder. This tool requires python3.4 

1) BAM to bedpe 

bamtobed tool from bedtools to convert bam to bedpe (bed file with paired end information)

bedtools bamtobed -bedpe -i <BAM>

https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html
  
2) sqlite database dump from bedpe

sql_dump.py to create the database from bedpe

3) Calculating the appropriate parameters for identifying topological boundaries

parametermcl.py to calculate the approriate bins and inflation value

4) Identifying the topological boundaries and classfying them based on their strength

mclclust.py to classify and identify boundaries
