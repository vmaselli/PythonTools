#!/usr/local/bin/perl58 -w
# Author:  Hye-Jung Elizabeth Chun
# Version: 1.0
# 2009-07-31

use strict;
use Data::Dumper;

use vars qw($opt_g $opt_i $opt_o $opt_n $opt_D);
use Getopt::Long qw(:config no_ignore_case);

GetOptions('g=s' => \$opt_g,
	   'i=s' => \$opt_i,
	   'o=s' => \$opt_o,
	   'n=s' => \$opt_n,
	   'D'   => \$opt_D
	   );

unless ($opt_g && $opt_i && $opt_o && $opt_n) {
  print "\n";
  print "create_WIG_file_from_SHRiMP_prob_outputlpl [OPTION]\n";
  print "\t-g  Original miRNA GFF file downloaded from the Sanger miRNA Registry\n";
  print "\t-i  Input SHRiMP prob file or Input directory where the SHRiMP prob output files are located\n";
  print "\t-D  Turn on if the intput is a directory\n";
  print "\t-o  Output gff file to be uploaded to the UCSC genome browser to see the mapping\n";
  print "\t-n  Name of the library and tray numbers (to be displayed as the name of the track)\n";
  print "\n";
  exit();
}


#Original miRNA GFF file example (format):
#4       .       miRNA   120442139       120442227       .       -       .       ACC="MI0000547"; ID="mmu-mir-30c-1";
#4       .       miRNA   120445211       120445302       .       -       .       ACC="MI0000259"; ID="mmu-mir-30e";
#4       .       miRNA   124408938       124409046       .       +       .       ACC="MI0004681"; ID="mmu-mir-697";
#4       .       miRNA   124421025       124421133       .       +       .       ACC="MI0004682"; ID="mmu-mir-698";

#SHRiMP probcalc file example (format):
##FORMAT: readname contigname strand contigstart contigend readstart readend readlength score editstring normodds pgenome pchance
#>1661_942_265_F3	mmu-mir-34a	+	20	54	1	35	35	2550	22x1A1x2x2x4x2	1.000000e+00	6.409997e-02	5.307948e-03
#>1659_1307_1823_F3	mmu-mir-24-2	+	61	92	1	32	35	2670	22x1x1T7	1.000000e+00	3.701575e-01	1.544394e-04
#>1662_2026_99_F3	mmu-mir-30c-1	-	40	73	1	34	35	2480	24x1x1A2G2x2	4.977171e-01	8.382310e-02	2.429946e-03

my %gff_table;
open(GFF, "<$opt_g") || die "Cannot open the file $opt_g\n";
while(<GFF>) {
  chomp;

  if($_ =~ /^\#/) {
    next;
  }

  my @lines = split "\t", $_;
  my $chr = $lines[0];
  my $chr_start = $lines[3];
  my $chr_end = $lines[4];
  my $id = $lines[8];

  if($id =~ /ID=\"(.+)\"/) {
    $id = $1;
  }

  $gff_table{$id}->{chr} = $chr;
  $gff_table{$id}->{chr_start} = $chr_start;
  $gff_table{$id}->{chr_end} = $chr_end;
}
close(GFF);



my %wig_table;
my $name = $opt_n;
open(OUTFILE, ">$opt_o");
print OUTFILE "track type=wiggle_0 name=\"${name}\" description=\"${name}\" visibility=full color=255,0,0 yLineMark=0 yLineOnOff=on priority=10\n";

if (defined $opt_D) {
  my @files = `ls $opt_i`;
  @files = map { chomp; my $tmp=$_ } @files;
  my @prob_files = grep { /.prob/ } @files;

  foreach my $file (@files) {
    $file = $opt_i.$file;
    open(PROB, "<$file");
    while(<PROB>) {
      chomp;

      if($_ =~ /^\#/) {
	next;
      }

      my @lines = split "\t", $_;
      my $id = $lines[1];
      my $contig_start = $lines[3];
      my $contig_end = $lines[4];

      my $miRNA_startPos = $gff_table{$id}->{chr_start};
      my $read_startPos = $miRNA_startPos + ($contig_start - 1);
      my $read_endPos = $read_startPos + ($contig_end - $contig_start);

      for (my $i=$read_startPos; $i<=$read_endPos; $i++) {
	$wig_table{$gff_table{$id}->{chr}}->{$i}++;
      }
    }
    close(PROB);
  }

  #Write the content of the wig hash table to the output file
  foreach my $chromosome (sort {$a<=>$b} keys %wig_table) {
    print OUTFILE "variableStep chrom=chr${chromosome} span=1\n";

    foreach my $base (sort {$a<=>$b} keys %{$wig_table{$chromosome}}) {
      print OUTFILE "$base\t", $wig_table{$chromosome}->{$base}, "\n";
    }
  }
}
else {
  open(PROB, "<$opt_i") || die "Cannot open the file $opt_i\n";

  while(<PROB>) {
    chomp;

    if($_ =~ /^\#/) {
      next;
    }

    my @lines = split "\t", $_;
    my $id = $lines[1];
    my $contig_start = $lines[3];
    my $contig_end = $lines[4];

    my $miRNA_startPos = $gff_table{$id}->{chr_start};
    my $read_startPos = $miRNA_startPos + ($contig_start - 1);
    my $read_endPos = $read_startPos + ($contig_end - $contig_start);

    for (my $i=$read_startPos; $i<=$read_endPos; $i++) {
      $wig_table{$gff_table{$id}->{chr}}->{$i}++;
    }
  }
  close(PROB);

  #Write the content of the wig has table to the output file
  foreach my $chromosome (sort {$a<=>$b} keys %wig_table) {
    print OUTFILE "variableStep chrom=chr${chromosome} span=1\n";

    foreach my $base (sort {$a<=>$b} keys %{$wig_table{$chromosome}}) {
      print OUTFILE "$base\t", $wig_table{$chromosome}->{$base}, "\n";
    }
  }
}

close(OUTFILE);

#Zip the file
`gzip $opt_o`;

print "Gzipping of the output file is done\n";
