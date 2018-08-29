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

perl (4 or higher)

python3

## Input file format
bam file; Sorted bam. (index file for bam (.bai) is required.)


microsatellte region file; List of microsatellte (tab-separated text or .gz file)
	Microsatellte region files are provided. But user defined microsatellte lies can be used. 

<chr> <start> <end> <repeat unit of microsatellte>
22      17283835        17283839        A
22      17283968        17283981        AT


reference genome file; Fasta file of reference genome (index file for samtools is required.)


parm.conf file; Parameter for microsatellte calling (optional)




## Output file format
<chr> <start> <end> <repeat unit of microsatellte> <number of reads with length of microsatellte> <genotype> <calling information (2nd major allele, number of reads, varinat allele frequency)>

## Usage
cd <path to MIV_VCALL>

perl RUN_MIV_CALL.pl -BAM <bam> -REF <reference.fas> -OUT <Output file name> -MS <microsatellte region file> -CONF <Configuration file (Optional)>

## Install
git clone @@@@@@@

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

##Contact

Akihiro Fujimoto - fujimoto@ddm.med.kyoto-u.ac.jp

https://sites.google.com/site/fujimotoakihironopeji/home/english