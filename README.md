# HiCMCL
High Performance Multi Resolution Markov Clustering for HiC data

All the required codes are in the bin folder. This tool requires python3.4 

1) BAM to bedpe

Here use bamtobed function from bedtools to convert bam to bedpe (bed file with paired end information)
bedtools bamtobed -bedpe -i <BAM>
  
2) sqlite database dump from bedpe

usage: sql_dump.py *args

----------SQL dump of bedpe file--------- [Date: 7th May 2018], [help: python sql_dump.py -h]

optional arguments:
  -h, --help            show this help message and exit
  
  -bedpe BEDPE, --bedpe BEDPE
                        comma seperated multiple bed file paths (single file can also be provided) containing paired interaction without header (Mandatory)
 
  -chrl CHRL, --chrl CHRL
                        chromosome length file containing two columns (chromosomenumber, length) without header (Mandatory)
 
  -cks CKS, --cks CKS   Chunk size for dumping to sql (for single iteration) (default=10000000)
  
  -mds MDS, --mds MDS   minimum intra chromosome distance (default=1000)
  
  -tag TAG, --tag TAG   comma seperated file tags (single tag can also be provided), the tags should match the -bedpe (Mandatory)
  
  -odir ODIR, --odir ODIR
                        outdir (Mandatory)
