Name
MIMcall
A software to indetify indels in microsatelte regions

Overview
1. Read extraction from cancer.bam and normal.bam from defined microsatellte regions 
2. Aanlysis the read information and calculation of likelihood using error rate matrix
3. Indetificaiton of indels in microsatellte regions

## Requirement
samtools
perl 
python3

## Usage
perl RUN_MIM_CALLER.pl -C_BAM <cancer.bam> -N_BAM <normal.bam> -REF <reference.fas> -OUT <Output file name> -MS <microsatellte regions>

## Install
git clone @@@@@@@

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## Author
Akihior Fujimoto