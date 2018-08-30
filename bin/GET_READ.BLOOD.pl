use strict;
use Getopt::Long;
use constant LF => "\n";
use constant {
    OPT_NO  => 0,
    OPT_YES => 1,
};

my $script_name = $0;
my @script_name_l = split("/", $script_name);
pop(@script_name_l);
my $SRC = join("/", @script_name_l);

my $mq_cutoff = 20;
my $len_cutoff1 = 100;
my $len_cutoff2 = 550;
my $flanking_len_cutoff = 3;
my $S_length_cutoff = 3;
my $q_score_cutoff = 10;
my $SW_range = 100;
#my $indel_S_num_cutoff = 2;
my ($match, $mismatch) = (1, 1);
my $SW_alignment = 0;
my $d = 1;
my $e = 1;
my ($BAM, $RMSK, $help, $bwa_ref) = ("", "", "", "", "");

GetOptions(
        "I=s" => \$BAM,
        "MS=s" => \$RMSK,
        "REF=s" => \$bwa_ref,
        "MQ=i" => \$mq_cutoff,
        "LL=i" => \$len_cutoff1,
        "ML=i" => \$len_cutoff2,
        "FL=i" => \$flanking_len_cutoff,
        "SL=i" => \$S_length_cutoff,
        "BQ=i" => \$q_score_cutoff,
#        "ISN=i" => \$indel_S_num_cutoff,
        "SW=i" => $SW_alignment,
        "GO=i" => \$d,
        "GE=i" => \$e,
        "h" => \$help
);

my $message = << '&EOT&';
Usage: perl <path>/GET_READ.BLOOD.pl [-I <input file (bam)>] [-MS <microsatellte location file>] [-REF <Reference file (fasta)>] [-h]
-I      Input bam file (Required)
-MS     MS location file (Required)
-REF    Reference genome (Index file for samtools is required.) (Required)
-h      Print this message
-MQ     Minimam mapping quality (20)
-LL     Minimum distance between reads (100)
-ML     Maximum distance between reads (550)
-FL     Minimum flanking length (3)
-SL     Minimum softclip length in a read (3)
-BQ     Minimum average base quality in softclip region (10)
-SW     Smith-Waterman alignment of read with softclip region (0)
-GO     Gap open penalty (1)
-GE     Gap extension penalty (1)
&EOT&

if($help == OPT_YES){
        print STDERR $message.LF;
        exit(0);
}

if(! $BAM or ! $RMSK){print STDERR $message.LF; exit(0);}

my %option_hash = ("I" => "/dev/stdin", "REF" => $bwa_ref, "MQ" => $mq_cutoff, "LL" => $len_cutoff1, "ML" => $len_cutoff2, "FL" => $flanking_len_cutoff, "SL" => $S_length_cutoff, "BQ" => $q_score_cutoff, "SW" => $SW_alignment, "GO" => $d, "GE" => $e);

my @option;
for my $name (keys %option_hash){
        my $tmp = "-"."$name"." "."$option_hash{$name}";
        push(@option, $tmp);
}
my $option_str = join(" ", @option);

open IN, "$RMSK" or die "Cannot open $RMSK";
while(<IN>){
	chomp;
	my @l = split("\t");
	my $chr = $l[0];
	my $start = $l[1];
	my $last = $l[2];

	my $MS_pos = "$chr".":".$start."-".$last;
	my $cmd = "samtools view -F 1024 -F 0x400 $BAM \"$MS_pos\"|perl $SRC/GPOS2RPOS.BLOOD.pl -CHR $chr -GP1 $start -GP2 $last $option_str";

	my %num;
	my $seq;
	foreach(`$cmd`){
		my @l = split("\t");
		$num{$l[5]}++;
		$seq = $l[6];
	}

	if(keys(%num) == 0){next;}
	print"$chr\t$start\t$last\t($l[3])n\t$seq\t";

	foreach(sort {$a <=> $b} keys %num){print"$_;$num{$_}\t";}
	print"\n";
}
