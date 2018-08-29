use strict;
use Getopt::Long;
use Data::Dumper;

our %error_rate;

my $blood_length_num_col = 6;
my $cancer_length_num_col = 5;
my $unit_col = 3;

my ($infile, $cancer_depth_cutoff, $blood_depth_cutoff, $blood_L_cutoff, $cancer_L_cutoff, $error_rate_file);

GetOptions(
	"I=s" => \$infile,
	"BD=i" => \$blood_depth_cutoff,
	"CD=i" => \$cancer_depth_cutoff,
	"BL=f" => \$blood_L_cutoff,
	"CL=f" => \$cancer_L_cutoff,
	"ER=s" => \$error_rate_file,
	"VAF=f" => \$min_VAF,
	"N=f" => \$min_num,

);

if($blood_depth_cutoff == undef){$blood_depth_cutoff = 10;}
if($blood_L_cutoff == undef){$blood_L_cutoff = -3;}
if($cancer_L_cutoff == undef){$cancer_L_cutoff = -3;}
if($min_VAF == undef){$min_VAF = 0.05;}
if($min_num == undef){$min_num = 2;}

my %error_rate_matrix = {};
open ER, "$error_rate_file" or  die "$error_rate_file !!";
my %header = {};
my $length_max = 0;
my $length_min = 0;
my $max_unit_length = 0;
while(<ER>){
	chomp;
	if($_ =~ /#/){
		my @l = split("\t");
		for(my $i = 2; $i < @l; $i++){
			$header{$i} = $l[$i];
		}
		$length_max = $l[-1];
		$length_min = $l[2];
		next;
	}

	my @l = split("\t");
	for(my $i = 2; $i < @l; $i++){
		$error_rate_matrix{$l[0]}->{$l[1]}->{$header{$i}} = $l[$i];
	}

	if($l[1] =~ /\//){
		my @unit = split("/", $l[1]);
		my $unit_length = length($unit[0]);
		if($max_unit_length < $unit_length){$max_unit_length = $unit_length}
	}
	else{
		my $unit_length = $l[1];
		$unit_length =~ s/bp//;
		if($max_unit_length < $unit_length){$max_unit_length = $unit_length}
	}
}

foreach my $range (keys %error_rate_matrix){
	foreach my $type (keys %{$error_rate_matrix{$range}}){
		my $total_error_rate = 0;
		foreach my $length (keys %{${$error_rate_matrix{$range}}{$type}}){
			$total_error_rate += $error_rate_matrix{$range}->{$type}->{$length};
		}
		if(exists($error_rate_matrix{$range}->{$type}->{0}) == 0){
			$error_rate_matrix{$range}->{$type}->{0} = 1 - $total_error_rate;
		}
	}
}

open IN, "$infile" or die "$infile !!";
while(<IN>){
	chomp;
#	print"$_\t";
	my @l = split("\t");
	my $pos = join("_", ($l[0], $l[1], $l[2]));

	my $unit = $l[$unit_col];
	$unit =~ s/\(//;
	$unit =~ s/\)n//;
	$unit = length($unit);

	my %error_rate = {};
	if($unit <= $max_unit_length){
		%error_rate = &make_error_rate(\%error_rate_matrix, $l[$unit_col], $l[2] - $l[1] + 1);
	}
	else{
		%error_rate = &make_error_rate(\%error_rate_matrix, "TTTTCT", $l[2] - $l[1] + 1);
	}

	my %blood_genotype2 = &get_blood_genotype2($l[$blood_length_num_col], $blood_L_cutoff, $blood_depth_cutoff, $unit, \%error_rate);

	if($blood_genotype2{"blood_genotype"} eq "LOW"){print"$_\tLOW\n"; next;}

	my $num_blood_allele = 0;
	my @tmp = split(",", $l[$blood_length_num_col]);
	$num_blood_allele = @tmp;

	my $cancer_total_depth = 0;
	foreach((split(",", $l[$cancer_length_num_col]))){
		my @tmp = split(";", $_);
		$cancer_total_depth += $tmp[1];
	}

	my $blood_total_depth = 0;
	foreach((split(",", $l[$blood_length_num_col]))){
		my @tmp = split(";", $_);
		$blood_total_depth += $tmp[1];
	}

	if($cancer_total_depth < $cancer_depth_cutoff){print"$_\tLOW\n"; next;}

	my @blood_genotype = split("/", $blood_genotype2{"blood_genotype"});
	@blood_genotype = sort { $a <=> $b } @blood_genotype;

	my $L_second_allele = "NA";
	if(exists($blood_genotype2{$blood_genotype[0]}) == 1 && $blood_genotype2{$blood_genotype[0]} != 1){$L_second_allele = $blood_genotype2{$blood_genotype[0]}}
	elsif(exists($blood_genotype2{$blood_genotype[1]}) == 1 && $blood_genotype2{$blood_genotype[1]} != 1){$L_second_allele = $blood_genotype2{$blood_genotype[1]}}

	my %blood_read_num;
	foreach((split(",", $l[$blood_length_num_col]))){
		my @tmp = split(";", $_);
		$blood_read_num{$tmp[0]} = $tmp[1];
	}

	my %cancer_allele = &get_cancer_specific_allele($blood_genotype2{"blood_genotype"}, $l[$cancer_length_num_col], $l[$blood_length_num_col], $cancer_L_cutoff, $unit);
	my $cancer_normal_allele_total;
	if(exists($cancer_allele{$blood_genotype[0]}) == 1 && exists($cancer_allele{$blood_genotype[1]}) == 1){
		$cancer_normal_allele_total = $cancer_allele{$blood_genotype[0]}->{read_num} + $cancer_allele{$blood_genotype[1]}->{read_num}
	}
	elsif(exists($cancer_allele{$blood_genotype[0]}) == 1){
		$cancer_normal_allele_total = $cancer_allele{$blood_genotype[0]}->{read_num};
	}
	elsif(exists($cancer_allele{$blood_genotype[1]}) == 1){
		$cancer_normal_allele_total = $cancer_allele{$blood_genotype[1]}->{read_num};
	}
	else{$cancer_normal_allele_total = 0;}

	if($blood_genotype[0] eq $blood_genotype[1]){
		my $L;
		foreach(keys %cancer_allele){
			if($cancer_allele{$_}->{germline} == 0){
				my $normal_allele_read_num = &set_normal_allele_read_num(\%cancer_allele, $blood_genotype[0], $cancer_normal_allele_total);					
				$L = &calculate_L_second($blood_genotype[0], $normal_allele_read_num, $_,  $cancer_allele{$_}->{read_num}, \%error_rate, $unit);
				if($L <= $cancer_L_cutoff){
					if(exists($blood_read_num{$_}) == 0){$blood_read_num{$_} = 0;}

					my $blood_L;
					if(exists($blood_genotype2{$_}) == 0){$blood_L = "NA"}
					else{$blood_L = $blood_genotype2{$_}}

					if($blood_L > $blood_L_cutoff){
						print"$_\tSignificant,$_;$cancer_allele{$_}->{read_num},$blood_read_num{$_},$L,$blood_L\n";
					}
				}
			}
		}
	}
	else{
		my $L;
                foreach(keys %cancer_allele){
			if($cancer_allele{$_}->{germline} == 0){
				my $L;
				my $normal_allele_read_num1 = &set_normal_allele_read_num(\%cancer_allele, $blood_genotype[0], $cancer_normal_allele_total);
				my $normal_allele_read_num2 = &set_normal_allele_read_num(\%cancer_allele, $blood_genotype[1], $cancer_normal_allele_total);
				for(my $i = 0; $i <= $cancer_allele{$_}->{read_num}; $i++){
					my ($cancer_allele1, $cancer_allele2) = ($i, $cancer_allele{$_}->{read_num} - $i);
					my $L1 = &calculate_L_second($blood_genotype[0], $normal_allele_read_num1, $_,  $cancer_allele1, \%error_rate, $unit);
					my $L2 = &calculate_L_second($blood_genotype[1], $normal_allele_read_num2, $_,  $cancer_allele2, \%error_rate, $unit);
					my $L1_2 = $L1 + $L2;
					if($L < -100){$L = -100}
					
					$L += 10**$L1_2;			
				}

				if($L == 0){$L = 0}
				else{				
					$L = &round(log($L)/log(10), 3);
				}
				if($L <= $cancer_L_cutoff){
					if(exists($blood_read_num{$_}) == 0){$blood_read_num{$_} = -1000;}

					my $blood_L;
					if(exists($blood_genotype2{$_}) == 0){$blood_L = "NA"}
					else{$blood_L = $blood_genotype2{$_}}
	
					if($blood_L > $blood_L_cutoff){
						print"$_\tSignificant,$_;$cancer_allele{$_}->{read_num},$blood_read_num{$_},$L,$blood_L\n";
					}
				}
			}
		}
	}

#	print"\n";
}

sub make_error_rate{
	my ($error_rate_matrix, $MS_type, $MS_range_length) = @_;
	my %error_rate_matrix = %{$error_rate_matrix};
	
	$MS_type =~ s/\(//;
	$MS_type =~ s/\)n//;
	my $MS_type_length = length($MS_type);
	
	my $found_range = 0;
	my $found_type = 0;
	my $range_in_matrix;
	my $MS_type_in_matrix;
	foreach my $range (keys %error_rate_matrix){
		my @range = split("_", $range);
		if($MS_range_length >= $range[0] && $MS_range_length <= $range[1]){
			$found_range = 1;
			$range_in_matrix = $range;
			$found_type = 0;
			foreach my $MS_type_in_matrix_tmp (keys %{$error_rate_matrix{$range}}){
				my @MS_type_in_matrix_tmp = ();
				if($MS_type_in_matrix_tmp =~ /\//){
					@MS_type_in_matrix_tmp = split("/", $MS_type_in_matrix_tmp);
					foreach(@MS_type_in_matrix_tmp){
						if($MS_type eq $_){
							$MS_type_in_matrix = $MS_type_in_matrix_tmp;
							$found_type = 1;
							last;
						}
					}
				}
				else{
					my $MS_type_tmp = length($MS_type);
					$MS_type_tmp .= "bp";
					if($found_type == 0 and $MS_type_tmp eq $MS_type_in_matrix_tmp){
						$MS_type_in_matrix = $MS_type_in_matrix_tmp;
						$found_type = 1;
					}
				}
			}
		}
	}

	if($found_range == 0 || $found_type == 0){
		die "Range or MS type did not found in the Error rate matrix!! $MS_type, $MS_range_length";
	}
	else{
		my %error_rate;
		foreach my $difference (keys %{$error_rate_matrix{$range_in_matrix}->{$MS_type_in_matrix}}){
			$error_rate{$MS_type_length}->{$difference} = $error_rate_matrix{$range_in_matrix}->{$MS_type_in_matrix}->{$difference};
		}
		return %error_rate;
	}
}

sub set_normal_allele_read_num{
	my ($cancer_allele, $normal_genotype, $cancer_normal_allele_total) = @_;
	my %cancer_allele = %{$cancer_allele};

	my $normal_alelle_read_num;
	if(exists($cancer_allele{$normal_genotype}->{read_num}) != 1 && $cancer_normal_allele_total > 0){$normal_alelle_read_num = int($cancer_normal_allele_total/2)}
	elsif(exists($cancer_allele{$normal_genotype}->{read_num}) != 1 && $cancer_normal_allele_total == 0){$normal_alelle_read_num = 1;}
	else{$normal_alelle_read_num = $cancer_allele{$normal_genotype}->{read_num};}

	return $normal_alelle_read_num;
}

sub get_cancer_specific_allele{
	my ($blood_genotype, $cancer_read_num, $blood_read_num, $cancer_L_cutoff, $unit) = @_;
	my @cancer_read_num = split(",", $cancer_read_num);
	my @blood_read_num = split(",", $blood_read_num);
	my @blood_genotype = split("/", $blood_genotype);

	my %cancer_allele;
	foreach(@cancer_read_num){
		my @cancer_tmp = split(";", $_);
		$cancer_allele{$cancer_tmp[0]}->{blood} = 0;
		$cancer_allele{$cancer_tmp[0]}->{germline} = 0;
		$cancer_allele{$cancer_tmp[0]}->{read_num} = $cancer_tmp[1];
	}

	foreach(@cancer_read_num){
		my @cancer_tmp = split(";", $_);
		foreach(@blood_read_num){
			my @blood_tmp = split(";", $_);
			if($cancer_tmp[0] == $blood_tmp[0]){$cancer_allele{$cancer_tmp[0]}->{blood} = 1;}
		}
		foreach(@blood_genotype){
			my @blood_tmp = split(";", $_);
			if($cancer_tmp[0] == $blood_tmp[0]){$cancer_allele{$cancer_tmp[0]}->{germline} = 1;}
		}
	}

	return %cancer_allele;
}

sub get_blood_genotype2{
	my ($blood_read_num, $L_cutoff, $depth_cutoff, $unit, $error_rate) = @_;
	my %major_second_alelle = &get_major_and_2nd($blood_read_num);
	my %error_rate = %{$error_rate};

	my %result;
	if($major_second_alelle{mn} + $major_second_alelle{sn} < $depth_cutoff){
                $result{"blood_genotype"} = "LOW";
		return %result;
        }
	elsif($major_second_alelle{sa} == 0){
		$result{"blood_genotype"} = "$major_second_alelle{ma}/$major_second_alelle{ma}";
		return %result;
        }
        elsif(($major_second_alelle{ma} - $major_second_alelle{sa})%$unit != 0){
	}
	
	my $L_second = &calculate_L_second($major_second_alelle{ma}, $major_second_alelle{mn}, $major_second_alelle{sa}, $major_second_alelle{sn}, \%error_rate, $unit);

	if($L_second > $L_cutoff){
		$result{"blood_genotype"} = "$major_second_alelle{ma}/$major_second_alelle{ma}";
		$result{$major_second_alelle{sa}} = $L_second;
		$result{$major_second_alelle{ma}} = 1;
	}
	else{
		$result{"blood_genotype"} = "$major_second_alelle{ma}/$major_second_alelle{sa}";
		$result{$major_second_alelle{sa}} = $L_second;
		$result{$major_second_alelle{ma}} = 1;
	}

	my @blood_genotype = split("/", $result{"blood_genotype"});
        @blood_genotype = sort { $a <=> $b } @blood_genotype;

	my @blood_read_num = split(",", $blood_read_num);	

	my %blood_genotype;
	foreach(@blood_read_num){
		my @tmp = split(";", $_);
		if($tmp[0] == $blood_genotype[0] || $tmp[0] == $blood_genotype[1]){
			$blood_genotype{$tmp[0]} = $tmp[1];
		}
	}

	foreach(@blood_read_num){
		my @tmp = split(";", $_);

		if($tmp[0] == $blood_genotype[0] || $tmp[0] == $blood_genotype[1]){next;}

		my $L;
		if(abs($tmp[0] - $blood_genotype[0]) < abs($tmp[0] - $blood_genotype[1])){
			$L = &calculate_L_second($blood_genotype[0], $blood_genotype{$blood_genotype[0]}, $tmp[0], $tmp[1], \%error_rate, $unit);
		}
		elsif(abs($tmp[0] - $blood_genotype[1]) < abs($tmp[0] - $blood_genotype[0])){
			$L = &calculate_L_second($blood_genotype[1], $blood_genotype{$blood_genotype[1]}, $tmp[0], $tmp[1], \%error_rate, $unit);
		}
		elsif(abs($tmp[0] - $blood_genotype[0]) == abs($tmp[0] - $blood_genotype[1])){
			if($tmp[0] < $blood_genotype[0]){
				$L = &calculate_L_second($blood_genotype[0], $blood_genotype{$blood_genotype[0]}, $tmp[0], $tmp[1], \%error_rate, $unit);
			}
			elsif($tmp[0] < $blood_genotype[1]){
				$L = &calculate_L_second($blood_genotype[1], $blood_genotype{$blood_genotype[1]}, $tmp[0], $tmp[1], \%error_rate, $unit);
			}
			elsif($tmp[0] > $blood_genotype[1]){
				$L = &calculate_L_second($blood_genotype[1], $blood_genotype{$blood_genotype[1]}, $tmp[0], $tmp[1], \%error_rate, $unit);
			}
			elsif($tmp[0] > $blood_genotype[0]){
				$L = &calculate_L_second($blood_genotype[0], $blood_genotype{$blood_genotype[0]}, $tmp[0], $tmp[1], \%error_rate, $unit);
			}
			else{warn "1 $tmp[0] $tmp[1] $blood_genotype[0] $blood_genotype[1] Something wrong ??";exit;}
		}
		else{warn "2 $tmp[0] $tmp[1] Something wrong ??";exit;}

		$result{$tmp[0]} = $L;
	}
	
	return %result;
}

sub round{
	my $A = shift;
	my $degit = shift;

	$A = int($A*(10**($degit - 1)) + 0.5)/(10**($degit - 1));
	return $A;
}

sub calculate_L_second{
	my ($ma, $mn, $sa, $sn, $error_rate, $unit) = @_;
	my %error_rate = %{$error_rate};

	my $A1 = &log_sum($mn + $sn);
	my $A2 = &log_sum($mn);
	my $A3 = &log_sum($sn);
	my $P_error_1 = $A1 - $A2 - $A3;
	my $difference;
	if(($sa - $ma)%$unit == 0){$difference = ($sa - $ma)/$unit;}
	else{
		if($sa - $ma > $unit){
			$difference = int(($sa - $ma)/$unit);
		}
		elsif($sa - $ma < $unit){
			$difference = 1;
		}
	}
	if($difference > $length_max){$difference = $length_max}
	elsif($difference < $length_min){$difference = $length_min}

	if($unit > $max_unit_length){$unit = $max_unit_length;}
	my $P_error_2 = $mn*(log(${$error_rate{$unit}}{0})/log(10));
	my $P_error_3 = $sn*(log(${$error_rate{$unit}}{$difference})/log(10));
	my $P_error = $P_error_1 + $P_error_2 + $P_error_3;
	$P_error = &round($P_error, 3);
	return $P_error;
}

sub get_major_and_2nd{
	my $allele = shift;
	my %allele;
	foreach((split(",", $allele))){
		my @tmp = split(";", $_);
		$allele{$tmp[0]} = $tmp[1];
	}

	my $order = 0;
	my $major_allele = "NA";
	my $major_num = 0;
	my $second_allele = "NA";
	my $second_num = 0;
	foreach(sort {$allele{$b} <=> $allele{$a}} keys %allele){
		$order++;
		if($order == 1){
			$major_allele = $_;
			$major_num = $allele{$_};
		}
		elsif($order == 2){
			$second_allele = $_;
			$second_num = $allele{$_};
		}
	}
	my %second_allele_dis;
	my $top_second_allele_dis = 0;
	my $top_second_allele_dis_num = 0;
	foreach(keys %allele){
		if($_ == $major_allele || $allele{$_} != $second_num){next;}
		$second_allele_dis{$_} = $_ - $major_allele;
		if(abs($second_allele_dis{$_}) > $top_second_allele_dis){$top_second_allele_dis = abs($second_allele_dis{$_});}
	}

	foreach(sort {$second_allele_dis{$a} <=> $second_allele_dis{$b}} keys %second_allele_dis){
		if($top_second_allele_dis == abs($second_allele_dis{$_})){
			$second_allele = $_;
			$second_num = $allele{$_};
		}
	}

	if($major_num == $second_num && $major_allele < $second_allele){
		my $tmp = $major_allele;
		$major_allele = $second_allele;
		$second_allele = $tmp;
	}
	my %major_second_alelle;
	%major_second_alelle = (ma => $major_allele, sa => $second_allele, mn => $major_num, sn => $second_num);
	return %major_second_alelle;
}

sub log_sum{
        my $A = shift;
        
        if($A == 0){return "log($A) cannot be calculated !!\n"; exit;}

        my $log_sum;
        for(my $i = 1; $i <= $A; $i++){$log_sum += log($i)/log(10);}

        return $log_sum;
}

