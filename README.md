Name
MIMcall
A software to indetify indels in microsatelte regions

Overview
1. Extract reads covering microsatellte regions from a bam file
2. Analyze repeat length 
3. Calculate likelihood using an error rate matrix, and indetificaiton of indels in microsatellte regions in cancer and matched nomal sample
4. Select cancer specifc miscoretellte

## Requirement
samtools (0.1.18 or higher)

perl (5 or higher)

python

## Input file format
**bam file**; Sorted bam (index file for bam (.bai) is required.)


**microsatellte region file**; List of microsatellte (tab-separated text or .gz file). Microsatellte region files are provided in MS_list directory. User defined microsatellte lists can be used. 

\<chr\> \<start\> \<end\> \<repeat unit of microsatellte\>  
22      17283835        17283839        A  
22      17283968        17283981        AT  


**parm.conf file**(optional) ; Parameter file for microsatellte calling (VIMcall/parm.conf is used. If you want to change paramters, please change this file (see below))


## Output file format
\<chr\> \<start\> \<end\> \<repeat unit of microsatellte\> \<Sequence of microsatellte> \<number of reads with length of microsatellte (length;number of reads)\> \<genotype\> \<calling information (2nd major allele, number of reads, varinat allele frequency)\>  
22      17282432        17282438        (A)n    AAAAAAA 7;44    7/7     -  
22      17282577        17282589        (A)n    AAAAAAAAAAAAA   13;24,14;3      13/14   minor_alelle=14;L=-3.03;Number=3;VAF=0.11 


## Usage
```
cd <path to MIVcall>
perl RUN_MIM_CALL.pl -C_BAM <cancer bam> -B_BAM <normal bam> -OUT <output directory> -MS <MS location file> -CONF <Config file (Optional)>
```

## Example
```
git clone https://github.com/afujimoto/MIMcall.git
cd MIMcall
perl RUN_MIM_CALL.pl  -C_BAM ./test/test.bam -B_BAM  ./test/test.2,bam -OUT ./test/test.out -MS./test/test_MS_list.txt
```


## Parameter setting in configuration file
We consider that the patemeter set of the provided configuration file is an apprppreate ones for 30x coverage WGS data. If you want to use different parameters, please change parm.conf file.

\##CANCER PRMS##  
mq_cutoff; Minimum mapping quality for reqd selection (20)  
len_cutoff1; Minimum distance between paired reads (100)  
len_cutoff2; Maxmum distance between paired reads (550)   
flanking_len_cutoff; Minimum flanking length (10)  
S_length_cutoff; Minimum softclip length (3)  
q_score_cutoff; Minimum average qility score of flanking region (10)  
indel_S_num_cutoff; Maximum number of indel and softclip. (2)
SW_alignment; Perform Smith-Waterman alignmnt (1; Yes, 0; No) (0)    
REF; Path to reference.fasta file (samtools index file is also required.)  
d; Gap open penalty (1)  
e; Gap extention penalty (1)  

\##NORMAL PRMS##
mq_cutoff; Minimum mapping quality for reqd selection (20)  
len_cutoff1; Minimum distance between paired reads (100)
len_cutoff2; Maxmum distance between paired reads (550) 
flanking_len_cutoff; Minimum flanking length (3)  
S_length_cutof; Minimum softclip length (3)
q_score_cutoff; Minimum average qility score of flanking region (10) 
SW_alignment; Perform Smith-Waterman alignmnt (1; Yes, 0; No) (0)
REF = "path to reference.fasta file (samtools index file is also required.)"
d; Gap open penalty (1)
e; Gap extention penalty (1)

\##MS CALL PRMS##		
BLOOD_MIN_DEPTH; Minimum depth for normal sample (15)  
CANCER_MIN_DEPTH; Minimum depth for cancer sample (15)  
BLOOD_L; Likelihood value for normal sample (-1)
CANCER_L; Likelihood value for cancer sample (-8)
ERROR_RATE_TABLE; Path of error arte matrix (VIMcall/Error_rate_matrix.txt)
VAF; Mimimum varinat allele frequency (0.05)    
NUM; Mimimum number of read (2)  

## Preformance
Performance of this tool is provided in Supplymanraty information of Fujimoto et al. (bioaxiv).

## Licence
GPL

## Contact

Akihiro Fujimoto - fujimoto@ddm.med.kyoto-u.ac.jp

https://sites.google.com/site/fujimotoakihironopeji/home/english

## Update