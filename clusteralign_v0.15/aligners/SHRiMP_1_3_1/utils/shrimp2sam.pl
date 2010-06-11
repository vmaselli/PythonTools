#!/bin/perl

# Contact: Nils Homer
# Version: 0.1.1

use strict;
use warnings;
use Getopt::Std;

&shrimp2sam;
exit;

sub shrimp2sam {
	my %opts = ();
	die("Usage: shrimp2sam.pl <shrimp pretty print .mapper>\n") if (@ARGV == 0 && -t STDIN);
	# core loop
	my (@s, $last, @staging, $k, $best_s, $subbest_s, $best_k);
	$last = '';
	while (<>) {
		if($_ =~ m/^>/) {
			my $l1 = $_;
			if(defined(my $l2 = <>) && defined(my $G = <>) && defined(my $l4 = <>) && defined(my $T = <>) && defined(my $R = <>)) {
				my ($name, $as) = &shrimp2sam_aux($_, $G, $T, $R, \@s); # read_name, number of mismatches
				if ($name eq $last) {
					# I do not know whether the multiple hits are ordered on the
					# number of mismatches. I assume they are not and so I have to
					# keep all these multiple hits in memory.
					@{$staging[$k]} = @s;
					if ($best_s > $as) {
						$subbest_s = $best_s;
						$best_s = $as;
						$best_k = $k;
					} elsif ($subbest_s > $as) {
						$subbest_s = $as;
					}
					++$k;
				} else {
					if ($last) {
						if ($best_s == $subbest_s) {
							$staging[$best_k][4] = 0;
						} elsif ($subbest_s - $best_s == 1) {
							$staging[$best_k][4] = 15 if ($staging[$best_k][4] > 15);
						}
						print join("\t", @{$staging[$best_k]}), "\n";
					}
					$k = 1; $best_s = $as; $subbest_s = -1; $best_k = 0;
					@{$staging[0]} = @s;
					$last = $name;
				}
			}
			else {
				die;
			}
		}
	}
#	print "best_k=$best_k\n";
	print join("\t", @{$staging[$best_k]}), "\n" if ($best_k >= 0);
}

sub shrimp2sam_aux {
	my ($line, $G, $T, $R, $s) = @_;
	chomp($line);
	chomp($G);
	chomp($T);
	chomp($R);
	my @t = split("\t", $line);
	my $ret;
	@$s = ();
	# read name
	$s->[0] = $ret = substr($t[0], 1);
	# initial flag (will be updated later)
	$s->[1] = 0;
	# read & quality
	$s->[9] = uc($T); $s->[9] =~ s/^T:\s+//g; $s->[9] =~ s/-//g;
	$s->[10] = "*";
	# cigar
	my $tmp = uc($t[9]); 
	$s->[5] = "";
	my $prev_m = 0;
	while(0 < length($tmp)) {
		if($tmp =~ m/^(\d+)/) { # ref match
			$prev_m += $1;
			$tmp = substr($tmp, length($1));
		}
		elsif($tmp =~ m/^([ACGTN])/) { # snp
			$prev_m++;
			$tmp = substr($tmp, 1);
		}
		elsif($tmp =~ m/^X/) { # crossover
			$tmp = substr($tmp, 1);
		}
		else { # indel 
			if(0 < $prev_m) {
				$s->[5] .= "$prev_m"."M"; $prev_m=0;
			}
			if($tmp =~ m/^\(([ACGTN]+)\)/) { # insertion
				$s->[5] .= "".length($1)."I";
				$tmp = substr($tmp, length($1) + 2);
			}
			elsif($tmp =~ m/^(\-+)/) { # deletion
				$s->[5] .= "".length($1)."D";
				$tmp = substr($tmp, length($1));
			}
			else {
				print "$tmp\n";
				die;
			}
		}
	}
	if(0 < $prev_m) {
		$s->[5] .= "$prev_m"."M"; $prev_m=0;
	}
	# coor
	$s->[2] = $t[1]; $s->[3] = $t[3];
	$s->[1] |= 0x10 if ($t[2] eq '-');
	# mapQ
	$s->[4] = 255;
	# ignore mate coordinate
	$s->[6] = '*'; $s->[7] = $s->[8] = 0;
	# aux
	my $as = $t[8];
	push(@$s, "AS:i:$as");
	return ($ret, $as);
}
