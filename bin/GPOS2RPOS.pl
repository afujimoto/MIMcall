use strict;
use Getopt::Long;
use constant LF => "\n";
use constant {
    OPT_NO  => 0,
    OPT_YES => 1,
};

my $mq_cutoff = 20;
my $len_cutoff1 = 100;
my $len_cutoff2 = 550;
my $flanking_len_cutoff = 10;
my $S_length_cutoff = 3;
my $q_score_cutoff = 10;
my $SW_range = 100;
my $indel_S_num_cutoff = 2;
my ($match, $mismatch) = (1, 1);
my $SW_alignment = 0;
my $d = 1;
my $e = 1;
my ($infile, $contig, $gpos1, $gpos1, $gpos2, $help, $bwa_ref) = ("", "", "", "", "", "", "");

GetOptions(
        "I=s" => \$infile,
	"CHR=s" => \$contig,
	"GP1=i" => \$gpos1,
	"GP2=i" => \$gpos2,
	"REF=s" => \$bwa_ref,
	"MQ=i" => \$mq_cutoff,
	"LL=i" => \$len_cutoff1,
	"ML=i" => \$len_cutoff2,
	"FL=i" => \$flanking_len_cutoff,
	"SL=i" => \$S_length_cutoff,
	"BQ=i" => \$q_score_cutoff,
	"ISN=i" => \$indel_S_num_cutoff,
	"SW=i" => $SW_alignment,
	"GO=i" => \$d,
	"GE=i" => \$e,
	"h" => \$help
);

my $message = << '&EOT&';
Usage: perl <path>/GPOS2RPOS.pl [-I <input file (bam)>] [-CHR <chr>] [-GP1 <MS start>] [-GP2 <MS end>] [-REF <Reference file (fasta)>] [-h]
-I	Input bam file (Required)
-CHR    Chr of MS (Required)
-GP1    Start of MS (Required)
-GP2	End of MS (Required)
-REF	Reference genome (Index file for samtools is required.) (Required)
-h	Print this message
-MQ	Minimam mapping quality (20)
-LL	Minimum distance between reads (100)
-ML	Maximum distance between reads (550)
-FL	Minimum flanking length (10)
-SL	Minimum softclip length in a read (3)
-BQ	Minimum average base quality in softclip region (10)
-ISN	Maximum number of indel and softclip in a read (2)
-SW	Smith-Waterman alignment of read with softclip region (0)
-GO	Gap open penalty (1)
-GE	Gap extension penalty (1)
&EOT&

if($help == OPT_YES){
	print STDERR $message.LF;
	exit(0);
}

open IN, "$infile" or die;
while(<IN>){
	chomp;
	my @l = split("\t");
	my $start = $l[3];
	my $cigar = $l[5];
	my $len = 0;
	my $total_len;
	my @type;
	my @len;

	my $filter = &read_filter2(\@l, $mq_cutoff, $len_cutoff1, $len_cutoff2);
	if($l[9] =~ /N/){
		$filter = 0; 
	}

	if($filter == 0){
		next;
	}

	my ($type, $len, $num_S, $total_len) = &cigar2len($cigar);
	@type = @{$type};
	@len = @{$len};

	my $last = $start + $total_len - 1;

	my $len_S = 0;
	if($num_S){
		for(my $i = 0; $i < @type; $i++){
			if($type[$i] eq "S"){$len_S += $len[$i];}	
		}
	}

	my $indel_S_num = 0;
        foreach(@type){
                if($_ eq "I" || $_ eq "D" || $_ eq "S"){$indel_S_num++;}
        }

	if($SW_alignment && (($num_S == 1 && $len_S > $S_length_cutoff) || ($indel_S_num >= $indel_S_num_cutoff))){
		my $S_q_value;
		if($num_S){
			$S_q_value = &S_q_value(\@type, \@len, $l[10]);
		}

		if(($num_S && $S_q_value > $q_score_cutoff) || ($indel_S_num >= $indel_S_num_cutoff)){
			my ($align_start_ref, $total_len);

			my ($ref_start, $ref_last) = &define_aligment_range($start, $last, \@type, \@len, $SW_range);
			my $ref_seq = &get_ref_seq($contig, $ref_start, $ref_last);
			my ($align_ref_seq1, $align_read_seq1, $align_start_ref1, $align_start_read1, $align_last_ref1, $align_last_read1) = &SW_align_gap_penalty($ref_seq, $l[9]);
			my ($cigar1, $mismatch1) = &align2cigar($align_ref_seq1, $align_read_seq1, $align_start_read1, $align_last_read1, length($l[9]));
			my @mismatch1 = @{$mismatch1};
			my ($type1, $len1, $num_S1, $total_len1) = &cigar2len($cigar1);
			my @type1 = @{$type1};
			my @len1 = @{$len1};

			my %type_num1 = &type_num(@type1);
			if($type_num1{S} > 0 || $type_num1{indel} > 1){
				my ($align_ref_seq2, $align_read_seq2, $align_start_ref2, $align_start_read2, $align_last_ref2, $align_last_read2) = &SW_align_gap_penalty(&rev_com($ref_seq), &rev_com($l[9]));
				$align_ref_seq2 = &rev_com($align_ref_seq2);
				$align_read_seq2 = &rev_com($align_read_seq2);
				$align_start_ref2 = length($ref_seq) - $align_last_ref2;
				$align_start_read2 = length($l[9]) - $align_last_read2;
				$align_last_ref2 = length($ref_seq) - $align_last_ref2 + &base_len($align_ref_seq2);
				$align_last_read2 = length($l[9]) - $align_last_read2 + &base_len($align_read_seq2);

				my ($cigar2, $mismatch2) = &align2cigar($align_ref_seq2, $align_read_seq2, $align_start_read2, $align_last_read2, length($l[9]));	
				my @mismatch2 = @{$mismatch2};
				my ($type2, $len2, $num_S2, $total_len2) = &cigar2len($cigar2);
				my @type2 = @{$type2};
				my @len2 = @{$len2};

				my %type_num2 = &type_num(@type2);
			
				if($type_num1{S} > $type_num2{S}){
					@type = @type2;
					@len = @len2;
					($align_start_ref, $total_len, $cigar) = ($align_start_ref2, $total_len2, $cigar2);
				}
				elsif($type_num2{S} >= $type_num1{S}){
					@type = @type1;
					@len = @len1;
					($align_start_ref, $total_len, $cigar) = ($align_start_ref1, $total_len1, $cigar1);
				}
				elsif($type_num1{indel} > $type_num2{indel}){
					@type = @type2;
					@len = @len2;
					($align_start_ref, $total_len, $cigar) = ($align_start_ref2, $total_len2, $cigar2);
				}
				elsif($type_num2{indel} >= $type_num1{indel}){
					@type = @type1;
					@len = @len1;
					($align_start_ref, $total_len, $cigar) = ($align_start_ref1, $total_len1, $cigar1);
				}
				else{
					if(@mismatch1 > @mismatch2){
						@type = @type2;
						@len = @len2;
						($align_start_ref, $total_len, $cigar) = ($align_start_ref2, $total_len2, $cigar2);
					}
					else{
						@type = @type1;
						@len = @len1;
						($align_start_ref, $total_len, $cigar) = ($align_start_ref1, $total_len1, $cigar1);
					}
				}
			}
			else{
				@type = @type1;
				@len = @len1;
				($align_start_ref, $total_len, $cigar) = ($align_start_ref1, $total_len1, $cigar1);
			}
		
			$start = $ref_start + $align_start_ref;
			$last = $ref_start + $align_start_ref + $total_len - 1;
		}
	}

	$indel_S_num = 0;
	foreach(@type){
		if($_ eq "I" || $_ eq "D" || $_ eq "S"){$indel_S_num++;}
	}
	if($indel_S_num >= $indel_S_num_cutoff){next;}

	my ($rpos1, $rpos2);
	if(($gpos1 - $flanking_len_cutoff >= $start && $gpos1 <= $last) && ($gpos2 >= $start && $gpos2 + $flanking_len_cutoff <= $last)){
		$rpos1 = &get_rpos(\@type, \@len, $start, $gpos1 - 1) + 1;
		$rpos2 = &get_rpos(\@type, \@len, $start, $gpos2 + 1) - 1;

		my $seq = substr($l[9], $rpos1, $rpos2 - $rpos1);
		if($rpos1 > 0 && $rpos2 < length($l[9]) - 1){
			my $ave_q1 = &ave_q_score(0, $rpos1 - 1, $l[10]);
			my $ave_q2 = &ave_q_score($rpos2 + 1, length($l[9]) - 1, $l[10]);
			if($ave_q1 < $q_score_cutoff || $ave_q2 < $q_score_cutoff){
				next;
			}
			
			my $seq_len = length($seq);

			print"$l[0]\t$l[1]\t$cigar\t$gpos1->$rpos1\t$gpos2->$rpos2\t$seq_len\t$seq\t$ave_q1\t$ave_q2\t$l[9]\t$_\n";
		}
		else{
		}
	}
}

sub rev_com{
	my $seq = shift;
	$seq =~ tr/[ATGC]/[TACG]/;
	$seq = reverse($seq);
	return $seq;
}

sub base_len{
	my $seq = shift;
	
	my $seq_len = 0;
	foreach((split("", $seq))){
		if($_ eq "-"){next;}
		$seq_len++;
	}

	return $seq_len;
}

sub type_num{
	my @type = @_;

	my %type_num = ("indel" => 0, "M" => 0, "S" => 0);
	foreach(@type){
		if($_ eq "D" || $_ eq "I"){$_ = "indel";}
		$type_num{$_}++;
	}
	
	return %type_num;
}

sub S_q_value{
	my ($type, $len, $q_val) = @_;
	my @type = @{$type};
	my @len = @{$len};

	my $ave_qval1;
	if($type[0] eq "S"){
		$ave_qval1 = &ave_q_score(0, $len[0] - 1, $q_val);
	}

	my $ave_qval2;
	if($type[-1] eq "S"){
		$ave_qval2 = &ave_q_score(length($q_val) - $len[-1], length($q_val) - 1, $q_val);
	}

	if($type[0] eq "S" && $type[-1] eq "S"){return(($ave_qval1 + $ave_qval2)/2);}
	elsif($type[0] eq "S"){return $ave_qval1;}
	elsif($type[-1] eq "S"){return $ave_qval2;}
	else{print"Somethig wrong in S_q_value !!"; exit;}
}

sub define_aligment_range{
	my ($start, $last, $type, $len, $SW_range) = @_;
	my @type = @{$type};
	my @len = @{$len};

	my ($ref_start, $ref_last);
	if($type[0] eq "S"){
		$ref_start = $start - $len[0] - $SW_range;
		$ref_last = $last;
	}
	elsif($type[-1] eq "S"){
		$ref_start = $start;
		$ref_last = $last + $len[-1] + $SW_range;
	}
	else{
		$ref_start = $start - $SW_range;
		$ref_last = $last + $SW_range;
	}
	
	return ($ref_start, $ref_last);
}

sub cigar2len{
	my $cigar = shift;

	my ($len, @type, @len);
	my $num_S = 0;
	my $total_len = 0;
	foreach($cigar =~ /(\d+)/g){
		my $num = $_;
		$len += length($num);
		my $type = substr($cigar, $len, 1);
		$len++;

		push(@type, $type);
		push(@len, $num);

		if($type eq "M"){$total_len += $num;}
		elsif($type eq "S"){$num_S++; next;}
		elsif($type eq "I"){next;}
		elsif($type eq "D"){$total_len += $num;}
	}

	return (\@type, \@len, $num_S, $total_len);
}

sub get_ref_seq{
	my ($contig, $ref_start, $ref_last) = @_;
	my $cmd = "samtools faidx $bwa_ref \"$contig:$ref_start-$ref_last\"";
	my $ref_seq;
	foreach(`$cmd`){
		chomp;
		if($_ =~ /^>/){next;}
		$ref_seq .= $_;
	}
	return $ref_seq;
}

sub SW_align_gap_penalty{
	my ($x, $y) = @_;
        my ($length_x, $length_y) = (length($x), length($y));

	my (@M, @X, @Y, @TRACE);
        for(my $i = 0; $i <= $length_x; $i++){
		$M[$i] = [(0) x ($length_y + 1)];
		$X[$i] = [(-100000) x ($length_y + 1)];
		$Y[$i] = [(-100000) x ($length_y + 1)];
	}

	for(my $i = 0; $i <= $length_x; $i++){
		for(my $j = 0; $j <= $length_y; $j++){
			for(my $k = 1; $k <= 3; $k++){#previous
				for(my $l = 1; $l <= 3; $l++){#current
					$TRACE[$i][$j][$k][$l] = 0;
				}
			}
		}
	}

	for(my $i=1; $i<=$length_x; $i++){
        	for(my $j=1; $j<=$length_y; $j++){
         		my $s;
			if(substr($x, $i - 1, 1) eq substr($y, $j - 1, 1)){$s = $match;}
			else{$s = (-1)*$mismatch;}

			my $val = &max(0, $M[$i-1][$j-1] + $s, $X[$i-1][$j-1] + $s, $Y[$i-1][$j-1] + $s);
			$M[$i][$j] = $val;
			if($val == $M[$i-1][$j-1] + $s){$TRACE[$i][$j][1][1] = 1;}
			if($val == $X[$i-1][$j-1] + $s){$TRACE[$i][$j][2][1] = 1;}
			if($val == $Y[$i-1][$j-1] + $s){$TRACE[$i][$j][3][1] = 1;}

			$val = &max($M[$i-1][$j] - $d, $X[$i-1][$j] - $e, $Y[$i-1][$j] - $d);
			$X[$i][$j] = $val;
			if($val == $M[$i-1][$j] - $d){$TRACE[$i][$j][1][2] = 1;}
			if($val == $X[$i-1][$j] - $e){$TRACE[$i][$j][2][2] = 1;}
			if($val == $Y[$i-1][$j] - $d){$TRACE[$i][$j][3][2] = 1;}

			$val = &max($M[$i][$j-1] - $d, $X[$i][$j-1] - $d, $Y[$i][$j-1] - $e);
			$Y[$i][$j] = $val;
			if($val == $M[$i][$j-1] - $d){$TRACE[$i][$j][1][3] = 1;}
			if($val == $X[$i][$j-1] - $d){$TRACE[$i][$j][2][3] = 1;}
			if($val == $Y[$i][$j-1] - $e){$TRACE[$i][$j][3][3] = 1;}
		}
	}

	my $max_score = 0;
	for(my $i=1; $i<=$length_x; $i++){
		for(my $j=1; $j<=$length_y; $j++){
			$max_score = &max($max_score, $M[$i][$j], $X[$i][$j], $Y[$i][$j]);
		}
	}

	my ($imax, $jmax, $lmax) = (0, 0, 0);
	for (my $i=1; $i<=$length_x; $i++){
		for (my $j=1; $j<=$length_x; $j++){
			if($max_score == $M[$i][$j]){
				$imax = $i;
				$jmax = $j;
				$lmax = 1;
			}
			if($max_score == $X[$i][$j]){
				$imax = $i;
				$jmax = $j;
				$lmax = 2;
			}
			if($max_score == $Y[$i][$j]){
				$imax = $i;
				$jmax = $j;
				$lmax = 3;
			}
		}
	}

	my ($Aligned_x, $Aligned_y, $align_start_x, $align_start_y) = &get_alignment($x, $y, \@TRACE, $imax, $jmax, $lmax);

	my ($align_last_x, $align_last_y) = ($imax, $jmax);

	return ($Aligned_x, $Aligned_y, $align_start_x, $align_start_y, $align_last_x, $align_last_y);
}

sub get_alignment{
	my ($x, $y, $TRACE, $i, $j, $l) = @_;
	my @TRACE = @{$TRACE};

	my ($Aligned_x, $Aligned_y) = ("", "");
	while($TRACE[$i][$j][1][$l] != 0 || $TRACE[$i][$j][2][$l] != 0 || $TRACE[$i][$j][3][$l] != 0){
		my $nextl;
		if($TRACE[$i][$j][3][$l]){$nextl = 3;}
		elsif($TRACE[$i][$j][2][$l]){$nextl = 2;}
		elsif($TRACE[$i][$j][1][$l]){$nextl = 1;}

		if($l == 1){
			$Aligned_x .= substr($x, $i - 1, 1);
			$Aligned_y .= substr($y, $j - 1, 1);
			$i--;
			$j--;
		}
		if($l == 2){
			$Aligned_x .= substr($x, $i - 1, 1);
			$Aligned_y .= "-";
			$i--;
		}
		if($l == 3){
			$Aligned_x .= "-";
			$Aligned_y .= substr($y, $j - 1, 1);
			$j--;
		}

		$l = $nextl;
	}
	
	$Aligned_x = reverse($Aligned_x);
	$Aligned_y = reverse($Aligned_y);

	my $align_start_x = $i;
	my $align_start_y = $j;
	return($Aligned_x, $Aligned_y, $align_start_x, $align_start_y);
}

sub max{
	my @matrix = @_;
	my $max = $matrix[0];

	foreach(@matrix){
		if($_ > $max){$max = $_;}
	}

	return $max;
}

sub SW_align{
	my ($seq1, $seq2) = @_;
	my $MATCH    =  1; # +1 for letters that match
	my $MISMATCH = -1; # -1 for letters that mismatch
	my $GAP      = -1; # -1 for any gap
	my $GAP_OPEN = -2;

	my @matrix;
	$matrix[0][0]{score} = 0;
	$matrix[0][0]{pointer} = "n";
	for(my $j = 1; $j <= length($seq1); $j++){
		$matrix[0][$j]{score} = 0;
		$matrix[0][$j]{pointer} = "n";
	}

	for(my $i = 1; $i <= length($seq2); $i++){
		$matrix[$i][0]{score} = 0;
		$matrix[$i][0]{pointer} = "n";
	}

	my $max_i = 0;
	my $max_j = 0;
	my $max_score = 0;

	for(my $i = 1; $i <= length($seq2); $i++){
		for(my $j = 1; $j <= length($seq1); $j++){
			my ($diagonal_score, $left_score, $up_score);
			#calculate match score
			my $letter1 = substr($seq1, $j - 1, 1);
			my $letter2 = substr($seq2, $i - 1, 1);
			if($letter1 eq $letter2){
				$diagonal_score = $matrix[$i - 1][$j - 1]{score} + $MATCH;
			}
			else{
				$diagonal_score = $matrix[$i - 1][$j - 1]{score} + $MISMATCH;
			}
			$up_score = $matrix[$i - 1][$j]{score} + $GAP;
			$left_score = $matrix[$i][$j - 1]{score} + $GAP;
			if($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0){
				$matrix[$i][$j]{score} = 0;
				$matrix[$i][$j]{pointer} = "n";
				next;
			}

			if($diagonal_score >= $up_score){
				if($diagonal_score >= $left_score){
					$matrix[$i][$j]{score} = $diagonal_score;
					$matrix[$i][$j]{pointer} = "d";
				}
				else{
					$matrix[$i][$j]{score} = $left_score;
					$matrix[$i][$j]{pointer} = "l";
				}
			}
			else{
				if($up_score >= $left_score){
					$matrix[$i][$j]{score} = $up_score;
					$matrix[$i][$j]{pointer} = "u";
				}
				else{
					$matrix[$i][$j]{score} = $left_score;
					$matrix[$i][$j]{pointer} = "l";
				}
			}

			if($matrix[$i][$j]{score} > $max_score){
				$max_i = $i;
				$max_j = $j;
				$max_score = $matrix[$i][$j]{score};
			}
		}
	}

	my $align1 = "";
	my $align2 = "";

	my $j = $max_j;
	my $i = $max_i;

	while(1){
		last if $matrix[$i][$j]{pointer} eq "n";
	
		if($matrix[$i][$j]{pointer} eq "d"){
			$align1 .= substr($seq1, $j - 1, 1);
			$align2 .= substr($seq2, $i - 1, 1);
			$i--;
			$j--;
		}
		elsif($matrix[$i][$j]{pointer} eq "l"){
			$align1 .= substr($seq1, $j - 1, 1);
			$align2 .= "-";
			$j--;
		}
		elsif($matrix[$i][$j]{pointer} eq "u"){
			$align1 .= "-";
			$align2 .= substr($seq2, $i - 1, 1);
			$i--;
		}
	}

	my $align_start1 = $j;
        my $align_start2 = $i;

	$align1 = reverse $align1;
        $align2 = reverse $align2;

	return ($align1, $align2, $align_start1, $align_start2);
}

sub align2cigar{
	my ($align1, $align2, $align1_start, $align1_last, $len_seq2) = @_;

	my $length2;
	if(length($align1) > length($align2)){
		$length2 = length($align1);
	}
	else{
		$length2 = length($align2);	
	}

	my $cigar;
        my $align_last = -1;
        my @mismatch;
        my $type;
        my %len;
        my $pre_type;
	my $align_len2 = 0;
	my $mismatch_pos = 0;

        for(my $j = 0; $j < $length2; $j++){
                my $letter1 = substr($align1, $j, 1);
                my $letter2 = substr($align2, $j, 1);

                if($letter1 =~ /[A-z]/){$align_last++;}
                if($letter2 =~ /[A-z]/){$align_len2++;}

                if($letter1 ne "-" && $letter2 ne "-"){
			$type = "M"; 
			$len{$type}++;
			$mismatch_pos++;
			if($letter1 ne $letter2){push(@mismatch, $mismatch_pos);}
		}
                elsif($letter1 eq "-" && $letter2 ne "-"){$type = "I"; $len{$type}++;}
                elsif($letter1 ne "-" && $letter2 eq "-"){$type = "D"; $len{$type}++; $mismatch_pos++;}

                if(defined($pre_type) == 1 && $pre_type ne $type){
                        $cigar .= "$len{$pre_type}"."$pre_type";
                        $len{$pre_type} = 0;
                }

		$pre_type = $type;
	}
	$cigar .= "$len{$pre_type}"."$pre_type";

	if(keys(%len) == 1){
		my $type_tmp = join("", keys(%len));
		$cigar = "$len{$type_tmp}"."$type_tmp";
	}

	if($align1_start != 0){$cigar = "$align1_start"."S"."$cigar";}
	if($len_seq2 > $align1_last){
		my $length_tmp = $len_seq2 - $align1_last;
		$cigar .= "$length_tmp"."S";
	}
	return ($cigar, \@mismatch);
}


sub get_mismatch_pos{
	my ($l, $cigar) = @_;
	my @l = @{$l};

	my %I_S_pos;

	my $len = 0;
	my $total_len;
	my @type;
	my @len;
	if($cigar =~ /I/ || $cigar =~ /S/){
		foreach($cigar =~ /[0-9]+[A-Z]+/g){
			my $num = $_;
			$num =~ s/[A-Z]+//;
			my $type = $_;
			$type =~ s/[0-9]+//;

			if($type eq "M"){$total_len += $num;}
			elsif($type eq "S"){$total_len += $num; $I_S_pos{$total_len} = $num;}
			elsif($type eq "I"){$total_len += $num; $I_S_pos{$total_len} = $num;}
			elsif($type eq "D"){next;}
		}
	}

	my $mismatch_info;
	foreach(@l){
		if($_ =~ /MD/){
				my @tmp = split(":", $_);
				$mismatch_info = $tmp[2];
		}
	}

	my @mismatch = split(/([^0-9]+)/, $mismatch_info);

	my $length;
	my @mismatch_pos;
	foreach(@mismatch){
		if($_ =~ /^\^/){next;}	
		$length += $_;
		if($_ =~ /[ATGC]/){
			$length++;
			push(@mismatch_pos, $length);
		}
	}

	if(keys(%I_S_pos) > 0 && @mismatch_pos > 0){
		foreach my $mismatch_pos (@mismatch_pos){
			foreach my $pos (keys %I_S_pos){
				if($pos < $mismatch_pos){$mismatch_pos += $I_S_pos{$pos};}
			}
		}
	}

	return @mismatch_pos;
}

sub read_filter{
	my ($l, $mq_cutoff, $len_cutoff1, $len_cutoff2) = @_;
	my @l = @{$l};

	my %filter;

	if($l[4] < $mq_cutoff){
		return 0;
	}

	foreach(@l){
		if($_ =~ /^XT/){ # alignment 
			my @tmp = split(":", $_);
			if($tmp[2] eq "R"){
				return 0;
			}
		}
		if($_ =~ /^Xt/){  #Not in original sam file
			my @tmp = split(":", $_);
			if($tmp[2] ne "U"){
				return 0;
			}
		}
		if($_ =~ /^MQ/){  #Not in original sam file
			my @tmp = split(":", $_);
			if($tmp[2] < $mq_cutoff){
				return 0;
			}
		}
		if($_ =~ /^XS/){ #Not in original sam file FR RF......
			my @tmp = split(":", $_);
			if($tmp[2] != 1){
				return 0;
			}
		}
		if($_ =~ /^Xl/){ #Not in original sam file
			my @tmp = split(":", $_);
			if($tmp[2] < $len_cutoff1 && $tmp[2] > $len_cutoff2){
				return 0;
			}
		}
	}
	
	return 1;
}

sub read_filter2{
	my ($l, $mq_cutoff, $len_cutoff1, $len_cutoff2) = @_;
	my @l = @{$l};
			
	if(abs($l[8]) < $len_cutoff1 || abs($l[8]) > $len_cutoff2){
		return 0;
	}

	my %filter;

	if($l[4] < $mq_cutoff){
		return 0;
	}

	if(! $l[1] & 0x2 ){
		return 0;
	}

	foreach(@l){
		if($_ =~ /^SA/){
			return 0;
		}
	}
	
	return 1;
}

sub ave_q_score{
	my ($start, $last, $q) = @_;
	
	my @q = split("", $q);
	foreach(@q){$_ = ord($_) - 33;}

	my $sum_q;
	for(my $i = $start; $i <= $last; $i++){$sum_q += $q[$i];}

	my $ave_q = $sum_q/($last - $start + 1);

	return $ave_q;
}


sub get_rpos{
	my ($type, $len, $start, $gpos) = @_;
	my @type = @{$type};
	my @len = @{$len};

	my @rpos;
	my $pre_gpos;
	my $pre_rpos;
	for(my $i = 0; $i < @type; $i++){
		if(! $pre_gpos){$pre_gpos = $start;}
		if($type[$i] eq "S"){
			for(my $j = 0; $j < $len[$i]; $j++){
				push(@{$rpos[0]}, "S");
				push(@{$rpos[1]}, $pre_gpos);
				push(@{$rpos[2]}, $pre_rpos + $j);
			}
			$pre_rpos = $pre_rpos + $len[$i];
		}
		elsif($type[$i] eq "M"){
			for(my $j = 0; $j < $len[$i]; $j++){
				push(@{$rpos[0]}, "M");
				push(@{$rpos[1]}, $pre_gpos + $j);
				push(@{$rpos[2]}, $pre_rpos + $j);
			}
			$pre_gpos = $pre_gpos + $len[$i];
			$pre_rpos = $pre_rpos + $len[$i];
		}
		elsif($type[$i] eq "I"){
			for(my $j = 0; $j < $len[$i]; $j++){
				push(@{$rpos[0]}, "I"); 
				push(@{$rpos[1]}, $pre_gpos);
				push(@{$rpos[2]}, $pre_rpos + $j);
			}
			$pre_rpos = $pre_rpos + $len[$i];
		}
		elsif($type[$i] eq "D"){
			for(my $j = 0; $j < $len[$i]; $j++){
				push(@{$rpos[0]}, "D"); 
				push(@{$rpos[1]}, $pre_gpos + $j);
				push(@{$rpos[2]}, $pre_rpos);
			}
			$pre_gpos = $pre_gpos + $len[$i];
		}
		else{die "Undefined cigar $type[$i] $start!!\n";}
	}

	my $rpos;
	for(my $k = 0; $k < @{$rpos[0]}; $k++){
		if($gpos == ${$rpos[1]}[$k]){
			$rpos = ${$rpos[2]}[$k];
		}
	}
	return $rpos;
}
