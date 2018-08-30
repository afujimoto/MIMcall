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

if(length($SRC) == 0){$SRC = "."}

my ($CANCER_BAM, $BLOOD_BAM, $OUTPUT_FILE, $RMSK, $config_file, $help, $bwa_ref) = ("", "", "","", "", "", "", "");

GetOptions(
        "C_BAM=s" => \$CANCER_BAM,
        "N_BAM=s" => \$BLOOD_BAM,
        "OUTPUT_F=s" => \$OUTPUT_FILE,
        "MS=s" => \$RMSK,
        "CONF=s" => \$config_file,
        "h" => \$help
);

my $message = << '&EOT&';
Usage: perl <path>/RUN_MIM_CALL.pl [-C_BAM <input file (cancer bam)>] [-B_BAM <input file (noarmal bam)>] [-OUT <output directory>] [-MS <MS location file>] [-CONF <Config file (Optional)>] [-h]
-C_BAM		Input cancer bam file (Required)
-N_BAM  	Input normal bam file (Required)
-MS		MS location file (Required)
-OUTPUT_F	Output file (Required)
-CONF		Config file (Optional)
-h      	Print this message
&EOT&

if($help == OPT_YES){
        print STDERR $message.LF;
        exit(0);
}

if(! $config_file){
	$config_file = "$SRC"."/"."parm.conf";
}

if( ! -f $CANCER_BAM){print"$CANCER_BAM Cancer bam file !!\n"; exit(0)}
if( ! -f $BLOOD_BAM){print"$BLOOD_BAM Normal bam file !!\n"; exit(0)}
if( ! $OUTPUT_FILE){print"$OUTPUT_FILE Outout file !!\n"; exit(0)}
if( ! -f $RMSK){print"$RMSK MS location file !!\n"; exit(0)}
if( ! -f $config_file){print"$config_file Config file !!\n"; exit(0)}

#####################GET PRMS###############################
my ($cancer, $normal, $call) = (0, 0, 0);
my %cancer_prms = {};
my %normal_prms = {};
my %call_prms = {};
open CONF, "$config_file" or die "Can not open $config_file !";
while(<CONF>){
	chomp;
	if($_ =~ /^#/ and $_ =~ /CANCER/){$cancer = 1; $normal = 0; $call = 0;}
	if($_ =~ /^#/ and $_ =~ /NORMAL/){$cancer = 0; $normal = 1; $call = 0;}
	if($_ =~ /^#/ and $_ =~ /CALL/){$cancer = 0; $normal = 0; $call = 1;}
	
	if ($_ !~ /=/){next;}

	my @l = split(/=/, $_);
	my $parm_name = $l[0];
	my $parm_val = $l[1];
	$parm_name =~ s/ //g;
	$parm_val =~ s/ //g;
#	print"parm_name /$parm_name/ parm_val /$parm_val/\n";
	
	if($cancer){$cancer_prms{$parm_name} = $parm_val}
	if($normal){$normal_prms{$parm_name} = $parm_val}
	if($call){$call_prms{$parm_name} = $parm_val}
}
############################################################

####################MAKE CMD################################
my %cancer_option = ("REF" => $cancer_prms{REF}, "MQ" => $cancer_prms{mq_cutoff}, "LL" => $cancer_prms{len_cutoff1}, "ML" => $cancer_prms{len_cutoff2}, "FL" => $cancer_prms{flanking_len_cutoff}, "SL" => $cancer_prms{S_length_cutoff}, "BQ" => $cancer_prms{q_score_cutoff}, "SW" => $cancer_prms{SW_alignment}, "GO" => $cancer_prms{d}, "GE" => $cancer_prms{e}, "MS" => $RMSK, "I" => $CANCER_BAM, "ISN" => $cancer_prms{indel_S_num_cutoff});

my %normal_option = ("REF" => $normal_prms{REF}, "MQ" => $normal_prms{mq_cutoff}, "LL" => $normal_prms{len_cutoff1}, "ML" => $normal_prms{len_cutoff2}, "FL" => $normal_prms{flanking_len_cutoff}, "SL" => $normal_prms{S_length_cutoff}, "BQ" => $normal_prms{q_score_cutoff}, "SW" => $normal_prms{SW_alignment}, "GO" => $normal_prms{d}, "GE" => $normal_prms{e}, "MS" => $RMSK, "I" => $BLOOD_BAM);

my $BC_merge_file = "$SRC"."/merged.txt";

#my %call_option = ("BD" => $call_prms{BLOOD_MIN_DEPTH}, "CD" => $call_prms{CANCER_MIN_DEPTH}, "BL" => $call_prms{BLOOD_L}, "CL" => $call_prms{CANCER_L}, "ER" => "$SRC/$call_prms{ERROR_RATE_TABLE}");
my %call_option = ("D" => $call_prms{MIN_DEPTH}, "L" => $call_prms{L}, "ER" => "$SRC/$call_prms{ERROR_RATE_TABLE}", "N" => "$call_prms{NUM}", "VAF" => "$call_prms{VAF}");

my @option = ();
for my $name (keys %cancer_option){
		if ($cancer_prms{SW_alignment} == 0 and $name eq "REF"){next;}
		if ($cancer_prms{SW_alignment} == 0 and $name eq "GO"){next;}
		if ($cancer_prms{SW_alignment} == 0 and $name eq "GE"){next;}
        my $tmp = "-"."$name"." "."$cancer_option{$name}";
        push(@option, $tmp);
}
my $cancer_option_str = join(" ", @option);

my @option = ();
for my $name (keys %normal_option){
		if ($normal_prms{SW_alignment} == 0 and $name eq "REF"){next;}
		if ($normal_prms{SW_alignment} == 0 and $name eq "GO"){next;}
		if ($normal_prms{SW_alignment} == 0 and $name eq "GE"){next;}
        my $tmp = "-"."$name"." "."$normal_option{$name}";
        push(@option, $tmp);
}
my $normal_option_str = join(" ", @option);

my @option = ();
for my $name (keys %call_option){
        my $tmp = "-"."$name"." "."$call_option{$name}";
        push(@option, $tmp);
}
my $call_option_str = join(" ", @option);

my $GET_READS_FROM_CANCER = "perl $SRC/bin/GET_READ.pl $cancer_option_str > $OUTPUT_FILE.CANCER";
my $GET_READS_FROM_NORMAL = "perl $SRC/bin/GET_READ.BLOOD.pl $normal_option_str > $OUTPUT_FILE.NORMAL";
my $MERGE = "python $SRC/bin/MERGE_MS.py $OUTPUT_FILE.CANCER $OUTPUT_FILE.NORMAL > $OUTPUT_FILE.merged";
my $CALL = "perl $SRC/bin/MIM_CALLER.pl -I $OUTPUT_FILE.merged $call_option_str > $OUTPUT_FILE";

############################################################

####################RUN#####################################
print"GET_READS_FROM_CANCER\n";
system("$GET_READS_FROM_CANCER");

print"GET_READS_FROM_NORMAL\n";
system("$GET_READS_FROM_NORMAL");

print"MERGE READS\n";
system("$MERGE");

print"MIM CALL\n";
system("$CALL");

#unlink $OUTPUT_FILE.CANCER;
#unlink $OUTPUT_FILE.NORMAL;
#unlink $OUTPUT_FILE.merged;

print"OUTPUT; $OUTPUT_FILE\n";
############################################################
