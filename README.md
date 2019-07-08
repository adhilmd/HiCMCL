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

3) Calculating the appropriate parameters for identifying topological boundaries

usage: parametermcl.py *args

-----------Parameter MCL-------- [Date: 7th May 2018], [help: python parametermcl.py -h]

optional arguments:
  -h, --help            show this help message and exit
  
  -findb FINDB, --findb FINDB
                        Comma seperated database file (Single file can also be provided) (Mandatory)
  
  -chrl CHRL, --chrl CHRL
                        Chromosome length file containing two columns (chromosomenumber, length) without header (Mandatory)
  
  -icore ICORE, --icore ICORE
                        Core resolution for clustering, in multiples (default=10)
  
  -nres NRES, --nres NRES
                        Base length for each random location (default=10000000)
  
  -cpath CPATH, --cpath CPATH
                        Main path were the codes are present for modules
  
  -norm NORM, --norm NORM
                        Normalization type 'median' or 'minmax' (default=minmax)
  
  -rand RAND, --rand RAND
                        Number of random location for the calculation (default=12)
  
  -th THREADS, --th THREADS
                        Number of threads (default=8)
  
  -tag TAG, --tag TAG   Comma seperated file tags for pickle file (single tag can also be provided), the tags should match the -bedpe
                        (Mandatory)
  
  -pref PREF, --pref PREF
                        Prefix for the output files (Mandatory)
  
  -odir ODIR, --odir ODIR
                        Outdir (Mandatory)

4) Identifying the topological boundaries and classfying them based on their strength

usage: mclclust.py *args


----------MCL clustering--------- [Date: 7th May 2018], [help: python mclclust.py -h]


optional arguments:
  -h, --help            show this help message and exit
  
  -findb FINDB, --findb FINDB
                        Comma seperated database file (Single file can also be provided) (Mandatory)
  
  -chrl CHRL, --chrl CHRL
                        Chromosome length file containing two columns (chromosomenumber, length) without header (Mandatory)
  
  -icore ICORE, --icore ICORE
                        Core resolution for clustering, in multiples (default=10)
  
  -norm NORM, --norm NORM
                        Normalization type 'median' or 'minmax' (default=minmax)
  
  -cpath CPATH, --cpath CPATH
                        Main path were the codes are present for modules
  
  -pfile PFILE, --pfile PFILE
                        parameter file containing resolution and inflation information from parametermcl script
  
  -th THREADS, --th THREADS
                        Number of threads (default=8)
  
  -tag TAG, --tag TAG   Comma seperated file tags for pickle file (single tag can also be provided), the tags should match the -bedpe
                        (Mandatory)
  
  -odir ODIR, --odir ODIR
                        Outdir (Mandatory)
