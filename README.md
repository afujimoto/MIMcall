# *MIMcall*

A software to indetify somatic indels in microsatellite regions

Overview
1. Extract reads covering microsatellte regions from cancer and normal bam files
2. Analyze repeat lengths 
3. Calculate likelihoods using an error rate matrix, and indetificaiton of indels in microsatellite regions in cancer and matched nomal sample

## Requirement
samtools (0.1.18 or higher)

perl (5 or higher)

python

## Input file format
**bam file**; Sorted bam (index file for bam (.bai) is required.)


**microsatellite region file**; List of microsatellite (tab-separated text or .gz file). Microsatellite region files are provided in MS_list directory. User defined microsatellite lists can be used. 

\<chr\> \<start\> \<end\> \<repeat unit of microsatellite\>  
22      17283835        17283839        A  
22      17283968        17283981        AT  


**parm.conf file**(optional) ; Parameter file for microsatellite calling (VIMcall/parm.conf is used by default. If you want to change parameters, please change this file (see below))


## Output file format
```
<chr> <start> <end> <repeat unit of microsatellite> <Sequence of microsatellite> <number of reads with length of microsatellite in cancer (length;number of reads)> <number of reads with length of microsatellite in normal tissue (length;number of reads)> <Calling result>  
1       159891271       159891275       (A)n    AAAAA   4;24,5;14,      4;17,5;12,      
1       42869711        42869720        (T)n    TTTTTTTTT       8;15,9;18,      8;1,9;14,10;1,          Mutant_allele=8;Number_in_cancer=15;Number_in_normal=1;L_of_cancer=-17.48;L_of_normal=-0.83;VAF=0.45    
```

Mutant_allele; Mutant allele  
Number_in_cancer; Number of reads of mutant allele in cancer  
Number_in_normal; Number of reads of mutant allele in normal tissue  
L_of_cancer; Likelihood of mutant allele in cancer  
L_of_normal; Likelihood of mutant allele in normal tissue  
VAF; Variant allele frequency in cancer  

## Usage
```
cd <path to MIMcall>
perl RUN_MIM_CALL.pl -C_BAM <cancer bam> -N_BAM <normal tissue bam> -OUT <output file> -MS <MS location file> -CONF <Config file (Optional)>
```

## Parameter setting in configuration file
We consider the parameter set in the provided configuration apprppreate for 30x coverage WGS data. If you would like to use different parameters, please make changes in the parm.config file. 

\##CANCER PRMS##   
mq_cutoff; Minimum mapping quality for read selection (20)   
len_cutoff1; Minimum distance between paired reads (100)   
len_cutoff2; Maximum distance between paired reads (550)    
flanking_len_cutoff; Minimum flanking length (10)   
S_length_cutoff; Minimum softclip length (3)   
q_score_cutoff; Minimum average quality score of flanking region (10)   
indel_S_num_cutoff; Maximum number of indel and softclip. (2)  

\##NORMAL PRMS##   
mq_cutoff; Minimum mapping quality for read selection (20)   
len_cutoff1; Minimum distance between paired reads (100)  
len_cutoff2; Maximum distance between paired reads (550)   
flanking_len_cutoff; Minimum flanking length (3)    
S_length_cutof; Minimum softclip length (3)  
q_score_cutoff; Minimum average quality score of flanking region (10)   

\##MS CALL PRMS##    
BLOOD_MIN_DEPTH; Minimum depth for normal sample (15)  
CANCER_MIN_DEPTH; Minimum depth for cancer sample (15)  
BLOOD_L; Likelihood value for normal sample (-1)  
CANCER_L; Likelihood value for cancer sample (-8)  
ERROR_RATE_TABLE; Path of error rate matrix (MIMcall/Error_rate_matrix.txt)  
BN; Maximum number of mutant reads in normal sample (1)  
CN; Minimum number of mutant reads in cancer (2)    
VAF; Minimum variant allele frequency in cancer (0.15)    

## Preformance
Performance of this tool is provided in Fujimoto et al. Genome Research (2020).
https://genome.cshlp.org/content/early/2020/03/16/gr.255026.119

## Licence
GPL

## Contact

Akihiro Fujimoto - afujimoto@m.u-tokyo.ac.jp

https://sites.google.com/site/fujimotoakihironopeji/home/english

## Update
