#!/usr/bin/perl
use strict;
use Getopt::Long;

my (
    $blastfile, $minidentity, $minlength, $maxmismatch,
    $maxgap,    $normnum,     $bymethods, $outfile
);
GetOptions(
    "bf=s"   => \$blastfile,
    "mi=s"   => \$minidentity,
    "ml=s"   => \$minlength,
    "mg=s"   => \$maxgap,
    "mm=s"   => \$maxmismatch,
    "norm=s" => \$normnum,
    "by=s"   => \$bymethods,
    "o=s"    => \$outfile
);

if(!defined $blastfile || !defined $minidentity || !defined $minlength || !defined $maxgap || !defined $maxmismatch || !defined $normnum || !defined $bymethods || !defined $outfile){
  die "Usage: GSMer.pl -bf <blastfile> -mi <minimum identity> -ml <minimum alignment length> -mg <maximum gap> -mm <maximum mismatch> -norm <normalization number> -by <normalization methods> -o <outfile>\n";
}

#print "$minidentity\t$minlength\t$maxgap\t$maxmismatch\t$normnum\t$bymethods\n";
my ( %sumhit, %sumblast, %hit, %blast );
open( BLAST, "$blastfile" ) || die "#1can not open $blastfile\n";
while (<BLAST>) {
    chomp;
    my @items = split( "\t", $_ );
    if (   $items[2] >= $minidentity
        && $items[3] >= $minlength
        && $items[4] <= $maxmismatch
        && $items[5] <= $maxgap )
    {
        $items[0] =~ /^(.*?)\_/;
        my $sample = $1;
        $items[1] =~ /^(\d+\_\d+)\-/;
        my $strainID = $1;
        $sumhit{$sample}{ $items[0] }                      = 1;
        $sumblast{$sample}{"$items[0]\_$items[1]"}         = 1;
        $hit{$strainID}{$sample}{ $items[0] }              = 1;
        $blast{$strainID}{$sample}{"$items[0]\_$items[1]"} = 1;
    }
}
close BLAST;
print keys %sumhit;

my %strain;
open(STRAIN,"/home/galaxy/galaxy-dist/tools/IEG/GSMer/strain.list")||die"#2 can not open strain.list\n";
while(<STRAIN>){
  chomp;
  my @items=split("\t",$_);
  $strain{$items[1]}=$items[0];
}
close STRAIN;

open( OUT, ">$outfile" ) || die "#3 can not open $outfile\n";
my @samples = keys %sumhit;
print OUT "strainID\tStrainName\t", join( "\t", @samples ), "\n";
foreach my $strainID ( keys %hit ) {
    print OUT "$strainID\t$strain{$strainID}";
    foreach my $sample (@samples) {
        my $num = scalar( keys %{ $hit{$strainID}{$sample} } )
          if $bymethods == 1;
        $num = scalar( keys %{ $blast{$strainID}{$sample} } )
          if $bymethods == 2;
          my $per = $num / scalar( keys %{ $sumhit{$sample} } ) * $normnum
          if $bymethods == 1;
        $per = $num / scalar( keys %{ $sumblast{$sample} } ) * $normnum if $bymethods == 2;
        print OUT "\t$per";
    }
    print OUT "\n";
}
close OUT;
          
