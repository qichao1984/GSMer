#
#===============================================================================
#
#         FILE: GSMer.pl
#
#  DESCRIPTION: This program identifies genome-specific markers from currently
#  sequenced micrbial genomes using k-mer based approaches.
#
#        FILES: GSMer.pl
#         BUGS: ---
#        NOTES: Exploring genome-specific markers for strain/species identification
#        in metagenomes
#       AUTHOR: Qichao Tu (qichao@ou.edu, philloid@gmail.com)
#  INSTITUTION: Institute for Environmental Genomics, University of Oklahoma
#      VERSION: 1.0
#      CREATED: 11/17/2013 10:28:38 AM
#     REVISION: ---
#===============================================================================
#!/usr/bin/perl
##INCLUDES
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Parallel::ForkManager;
use String::Random;
require 'Configuration.pm';

######################
##VARIABLES
######################
my ( $exmode, $list, $infile );

GetOptions(
  "m=s"  => \$exmode,
  "l=s"  => \$list,
  "i=s"  => \$infile,
  "f1=s" => \$inf1file
);

######################
##PARAMETERS
######################
#programs
my $formatdb = $Configuration::formatdb || die "Please specify formatdb path\n";
my $meryl    = $Configuration::meryl    || die "Please specify meryl path\n";
my $mapMers  = $Configuration::mapMers  || die "Please specify mapMers path\n";
my $blastall = $Configuration::blastall || die "Please specify blastall path\n";

##GSM parameters
#taxonomic level, strain or species
my $tax = $Configuration::tax || die "Please specify taxonomic level\n";

#k-mer size for stretch filtering
my $k = $Configuration::k || die "Please specify kmer size\n";

#GSM length
my $probe_length = $Configuration::probe_length || die "Please specify probe length\n";

#continuous stretch cutoff
my $stretch_cutoff = $Configuration::stretch_cutoff || die "Please specify stretch cutoff\n";

#identity cutoff
my $identity_cutoff = $Configuration::identity_cutoff || die "Please specify identity cutoff\n";

##output directories
#root directory for output files
my $root_directory = $Configuration::root_directory
  || die "Please specify root directory\n";

#directory for storing genbank files for strain or species
my $gbkdir = $root_directory . "/" . $Configuration::gbkdir
  || die "Please specify gbk file directory\n";
mkdir $gbkdir unless ( -e $gbkdir );    #create gbk fold if not exist

#directory for storing blast database
my $blastdbdir = $root_directory . "/" . $Configuration::blastdbdir
  || die "Please specify gbk file directory\n";
mkdir $blastdbdir unless ( -e $blastdbdir );

#directory for stroing k-mer databases
my $merdbdir = $root_directory . "/" . $Configuration::merdbdir
  || die "Please specify merdb directory\n";
mkdir $merdbdir unless ( -e $merdbdir );

#directory for output GSMs
my $probedir = $root_directory . "/" . $Configuration::probedir
  || die "Please specify probe directory\n";
mkdir $probedir unless ( -e $probedir );

##multi-threads
#number of threads be used in the program
my $threads = $Configuration::threads || die "Please specify threads number\n";

######################
##BEGIN
######################
BEGIN:{
  &print_start();
}
######################
##MAIN
######################
MAIN: {
  ##help and check input options
  &printhelp() if !$exmode;
  &printhelp() if $exmode eq "help";
  if ( $exmode ne "splitgbk"
    && $exmode ne "readgbk"
    && $exmode ne "makeblastdb"
    && $exmode ne "makekmerdb"
    && $exmode ne "getgsm"
    && $exmode ne "mapgsm"
    && $exmode ne "blastgsm"
    && $exmode ne "checkspecificity" )
  {
    die "unrecognized execution mode option: $exmode\n";
  }
  ##split combined gbk file into strain/species level
  if ( $exmode eq "splitgbk" ) {
    if ( !$infile ) {
      die "Missing input genbank file!\n";
    }
    print "Spliting $infile into strain/species level gbk files...\n";
    %genomelist = &SplitGBKByTax( $infile, $gbkdir, $tax );
    print "A total of ",scalar(keys %genomelist)," gbk files were generated!\n";
    open( LIST, ">strain.list" ) || die "# can not write strain.list";
    foreach my $strain ( keys %genomelist ) {
      print LIST "$strain\t$genomelist{$strain}\n";
    }
    close LIST;
  }

  ##read individual gbk file from gbk directory and generate strain.list file
  if ( $exmode eq "readgbk" ) {
    print "Reading genbank files from $gbkdir...\n";
    %genomelist = &ReadTaxFromGBK( $gbkdir, $tax );
    print "A total of ",scalar(keys %genomelist)," gbk files were read!\n";
    open( LIST, ">strain.list" ) || die "# can not write strain.list";
    foreach my $strain ( keys %genomelist ) {
      print LIST "$strain\t$genomelist{$strain}\n";
    }
    close LIST;
  }

  ##make blast database from all microbial genomes and/or alien genomes (e.g. human genome)
  if ( $exmode eq "makeblastdb" ) {
    my $blastdbfile = "$blastdbdir/f1f2.all.fa";
    print "Generating blast database from microbial genomes and/or alien genomes...\n";
    open( FA, ">$blastdbfile" ) || die "# can not write $blastdbfile\n";
    foreach my $gbk ( glob("$gbkdir/*gbk") ) {
      $gbk =~ /$gbkdir\/(.*?).gbk/;
      my $i      = $1;
      my $infile = "$gbkdir/${i}.gbk";
      my $seqIO  = Bio::SeqIO->new( -file => "$infile", -format => "GenBank" );
      while ( my $seqobj = $seqIO->next_seq ) {
        my $seq = $seqobj->seq;
        my $id  = $seqobj->id;
        print FA ">$i\-$id\n$seq\n";
      }
    }
    if ( $inf1file ne "" ) {
      my $f1seqIO = Bio::SeqIO->new( -file => "$inf1file", -format => "fasta" );
      while ( my $seqobj = $f1seqIO->next_seq ) {
        my $seq = $seqobj->seq;
        my $id  = $seqobj->id;
        print FA ">$id\n$seq\n";
      }
    }
    close FA;
    system("$formatdb -i $blastdbfile -p F");
    print "Done formatting blast database $blastdbfile!\n";
  }

  ##make k-mer database from all microbial genomes and/or alien genomes (e.g. human genome)
  if ( $exmode eq "makekmerdb" ) {
    print "Generating k-mer database from microbial genomes...\n";
    my $pm = new Parallel::ForkManager($threads);
    my %gbk;
    foreach my $gbk ( glob("$gbkdir/*gbk") ) {
      $gbk =~ /$gbkdir\/(.*?).gbk/;
      my $ij = $1;
      if ( $tax == 0 ) {
        push( @{ $gbk{$ij} }, $gbk );
      }
      if ( $tax == 1 ) {
        $ij =~ /(\d+)\_(\d+)/;
        my $i = $1;
        push( @{ $gbk{$i} }, $gbk );
      }
    }

    foreach my $ij ( keys %gbk ) {
      my $pid        = $pm->start and next;
      my @gbk        = @{ $gbk{$ij} };
      my $outfile    = "$gbkdir/${ij}.${k}mer.fa";
      my $kmerdbfile = "$gbkdir/${ij}.${k}mer";
      ##generate k-mer fasta files for each strain
      if ( !-e $outfile ) {
        &GenerateKmers( \@gbk, $outfile, "GenBank", $k, $ij );
        &GenerateKmerDB( $outfile, $kmerdbfile, $meryl, 10, $threads, $k );
      }
      $pm->finish;
    }
    $pm->wait_all_children;
    
    ##merge kmer db
    print "Merging strain/species level k-mer tables...\n";
    my $f1kmerdb     = "$merdbdir/f1all.k${k}";
    my $f2kmerdb     = "$merdbdir/f2all.k${k}";
    my $f2kmerdbn2fa = "$merdbdir/f2n2.k${k}.fa";
    my $f2kmerdbn2   = "$merdbdir/f2n2.k${k}";
    my $f1f2n2kmerdb = "$merdbdir/f2n2_f1.k${k}";

    my @kmerdbfiles;
    foreach my $file ( glob("$gbkdir/*${k}mer.mcdat") ) {
      $file =~ s/.mcdat//;
      push( @kmerdbfiles, $file );
    }
    if ( !-e "$merdbdir/$f2kmerdb\.mcdat" ) {
      if ( $#kmerdbfiles > 99 ) {
        my $t = int( $#kmerdbfiles + 1 ) / 100;
        my @tmpfiles;
        for ( my $i = 0; $i <= $t; $i++ ) {
          my $tmpfile = "tmp$i";
          push( @tmpfiles, $tmpfile );
          my $meryladdcmd = "$meryl -M add";
          for ( my $j = $i * 100; $j <= ( $i + 1 ) * 100 - 1; $j++ ) {
            my $merfile = $kmerdbfiles[$j];
            $meryladdcmd .= " -s $merfile" if $merfile;
          }
          $meryladdcmd .= " -o $tmpfile";
          system("$meryladdcmd");
        }

        my $meryladdcmd = "$meryl -M add";
        foreach my $tmpfile (@tmpfiles) {
          $meryladdcmd .= " -s $tmpfile";
        }
        $meryladdcmd .= " -o $f2kmerdb";
        system("$meryladdcmd");

        foreach my $tmpfile (@tmpfiles) {
          unlink "$tmpfile.*";
        }
      }
      else {
        my $meryladdcmd = "$meryl -M add";
        foreach my $file (@kmerdbfiles) {
          $meryladdcmd .= " -s $file";
        }
        $meryladdcmd .= " -o $f2kmerdb";
        system("$meryladdcmd");
      }
    }

    ##dump k-mers with frequency>=2, i.e. showing up in two or more strains
    if ( !-e $f2kmerdbn2fa ) {
      print "Extracting k-mers showing up in >=2 microbial genomes...\n";
      system("$meryl -Dt -n 2 -s $f2kmerdb > $f2kmerdbn2fa");
      &GenerateKmerDB( $f2kmerdbn2fa, $f2kmerdbn2, $meryl, 10, $threads, $k );
    }

    ##generate k-mer database for f1file that all k-mers were kept, combine with f2 kmer databse with frequency cutoffs
    if ( $inf1file ne "" ) {
      if ( !-e "$f1f2n2kmerdb.mcdat" ) {
        print "Adding alien k-mers into the k-mer database...\n";
        &GenerateKmerDB( $inf1file, $f1kmerdb, $meryl, 10, $threads, $k );
        system("$meryl -M add -s $f2kmerdbn2 -s $f1kmerdb -o $f1f2n2kmerdb");
      }
    }
  }
  ##generate all candidate GSMs for interested microbial genomes
  if ( $exmode eq "getgsm" ) {
    print "Generating candidate $probe_length bp genome-specific markers...\n";
    my %genomelist;
    open( LIST, "strain.list" ) || die "# can not open strain.list";
    while (<LIST>) {
      chomp;
      my @items = split( "\t", $_ );
      $genomelist{ $items[0] } = $items[1];
    }
    close LIST;

    my %list;
    if ($list) {
      open( LIST, "$list" ) || die "# can not open $list\n";
      while (<LIST>) {
        chomp;
        $list{$_} = 1 if $genomelist{$_};
      }
      close LIST;
    }

    my %glist = %list if $list;
    %glist = %genomelist if !$list;
    my $pm = new Parallel::ForkManager($threads);
    foreach my $strain ( keys %glist ) {
      my $pid          = $pm->start and next;
      my $i            = $genomelist{$strain};
      my $infile       = "$gbkdir/${i}.gbk";
      my $probemerfile = "$probedir/${i}.${probe_length}mer.fa";
      &GenerateGSMs( $infile, $probemerfile, "GenBank", $probe_length, $i );
      $pm->finish;
    }
    $pm->wait_all_children;
  }
  ##map candidate GSMs to k-mer database for continuous stretch filtering
  if ( $exmode eq "mapgsm" ) {
    my $f1f2n2kmerdb;
    print "Mapping GSMs to k-mer database for continuous stretch filtering...\n";
    #select which k-mer database should be used
    if ( -e "$merdbdir/f2n2_f1.k${k}.mcdat" ) {
      $f1f2n2kmerdb = "$merdbdir/f2n2_f1.k${k}";
    }
    else {
      $f1f2n2kmerdb = "$merdbdir/f2n2.k${k}";
    }

    #get genome list from strain.list file
    my %genomelist;
    open( LIST, "strain.list" ) || die "# can not open strain.list";
    while (<LIST>) {
      chomp;
      my @items = split( "\t", $_ );
      $genomelist{ $items[0] } = $items[1];
    }
    close LIST;

    #read interested microbial strains
    #if no list was supplied, all genomes in strain.list will be used
    my %list;
    if ($list) {
      open( LIST, "$list" ) || die "# can not open $list\n";
      while (<LIST>) {
        chomp;
        $list{$_} = 1 if $genomelist{$_};
      }
      close LIST;
    }
    else {
      %list = %genomelist;
    }

    #map candidate GSMs to k-mer database for continuous strech filtering
    #multi-threads is enabled
    my $pm = new Parallel::ForkManager($threads);
    foreach my $strain ( keys %list ) {
      my $pid           = $pm->start and next;
      my $strainID      = $genomelist{$strain};
      my $probemerfile  = "$probedir/${strainID}.${probe_length}mer.fa";
      my $probemermap   = $probemerfile;
      my $probemermapfa = $probemerfile;
      $probemermap   =~ s/fa$/map${k}/;
      $probemermapfa =~ s/fa$/map${k}.fa/;
      my $size = -s $probemermap;
      &MapAndFilter( $f1f2n2kmerdb, $probemerfile, $probemermap, $probemermapfa, $k )
        if ( $size == 0 );
      $pm->finish;
    }
    $pm->wait_all_children;
  }
  ##run blast search unmapped GSMs for similarity based filtering
  if ( $exmode eq "blastgsm" ) {
    my $blastdbfile = "$blastdbdir/f1f2.all.fa";
    print "Blast searching GSMs againt blast database...\n";
    my %genomelist;
    open( LIST, "strain.list" ) || die "# can not open strain.list";
    while (<LIST>) {
      chomp;
      my @items = split( "\t", $_ );
      $genomelist{ $items[0] } = $items[1];
    }
    close LIST;

    my %list;
    if ($list) {
      open( LIST, "$list" ) || die "# can not open $list\n";
      while (<LIST>) {
        chomp;
        $list{$_} = 1 if $genomelist{$_};
      }
      close LIST;
    }
    else {
      %list = %genomelist;
    }
    my $pm = new Parallel::ForkManager($threads);
    foreach my $strain ( keys %list ) {
      my $pid               = $pm->start and next;
      my $strainID          = $genomelist{$strain};
      my $probemermapfa     = "$probedir/${strainID}.${probe_length}mer.map${k}.fa";
      my $probemermapsblast = $probemermapfa;
      $probemermapsblast =~ s/fa$/sblast/;
      my $line = `grep -c ">" $probemermapfa`;
      if ( $line >= 100 ) {

        #megablast is used here
        #change the parameters if blastall or other blast like programs are used
        system(
          "$blastall -d $blastdbfile -i $probemermapfa -o $probemermapsblast -m 8 -F F -p 90 -W 15"
        );
      }
      $pm->finish;
    }
    $pm->wait_all_children;
  }
  ##check GSM specificity
  if ( $exmode eq "checkspecificity" ) {
    print "Cheking GSM specificity based on blast results...\n";
    my $blastdbfile = "$blastdbdir/f1f2.all.fa";
    my %genomelist;
    open( LIST, "strain.list" ) || die "# can not open strain.list";
    while (<LIST>) {
      chomp;
      my @items = split( "\t", $_ );
      $genomelist{ $items[0] } = $items[1];
    }
    close LIST;

    my %list;
    if ($list) {
      open( LIST, "$list" ) || die "# can not open $list\n";
      while (<LIST>) {
        chomp;
        $list{$_} = 1 if $genomelist{$_};
      }
      close LIST;
    }
    else {
      %list = %genomelist;
    }

    my $pm = new Parallel::ForkManager($threads);
    foreach my $strain ( keys %list ) {
      my $pid               = $pm->start and next;
      my $strainID          = $genomelist{$strain};
      my $probemermapfa     = "$probedir/${strainID}.${probe_length}mer.map${k}.fa";
      my $probemermapsblast = $probemermapfa;
      $probemermapsblast =~ s/fa$/sblast/;
      my $tmpprobeoutfile = $probemermapfa;
      $tmpprobeoutfile =~ s/fa$/out/;
      &CheckProbeSpecificity( $probemermapfa, $probemermapsblast, $blastdbfile,
        $tmpprobeoutfile )
        if ( -e $probemermapsblast && !-e $tmpprobeoutfile );
      $pm->finish;
    }
    $pm->wait_all_children;
  }
}

###################################################################################
sub SplitGBKByTax() {
###################################################################################
  my ( $gbkfile, $outdir, $tax ) = @_;
  my %list;
  my %specieslist;
  my %strainlist;
  my $seqIO = Bio::SeqIO->new( -file => "$gbkfile", -format => "GenBank" );
  while ( my $seqobj = $seqIO->next_seq ) {
    my @organism    = $seqobj->species->classification;
    my $strain      = $organism[$tax];
    my @strainitems = split( " ", $strain );
    my $species     = "$strainitems[0] $strainitems[1]";
    my $outgbk;
    if ( !$specieslist{$species} ) {
      $specieslist{$species} = scalar( keys %specieslist ) + 1;
    }
    if ( !$strainlist{$species}{$strain} ) {
      my $i = $specieslist{$species};
      my $j = scalar( keys %{ $strainlist{$species} } ) + 1;
      $outgbk                        = "${i}_${j}.gbk";
      $list{$strain}                 = "${i}_${j}";
      $strainlist{$species}{$strain} = 1;
    }
    else {
      my $ij = $list{$strain};
      $outgbk = "${ij}.gbk";
    }
    my $out = Bio::SeqIO->new( -file => ">>$outdir/$outgbk", -format => "GenBank" );
    $out->write_seq($seqobj);
  }
  return %list;
}

###################################################################################
sub ReadTaxFromGBK() {
###################################################################################
  my ( $dir, $tax ) = @_;
  my %genomelist;
  my $pm = new Parallel::ForkManager( $threads * 4 );
  $pm->run_on_finish(
    sub {
      my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $tmplist ) = @_;
      if ( defined($tmplist) ) {
        %genomelist = ( %genomelist, %$tmplist );
      }
      else {
        print qq|No message received from child process $pid!\n|;
      }
    }
  );
  foreach my $file ( glob("$dir/*gbk") ) {
    my $pid = $pm->start and next;
    $file =~ /(.*?)\.gbk/;
    my $i        = $1;
    my $seqIO    = Bio::SeqIO->new( -file => "$file", -format => "GenBank" );
    my $seqobj   = $seqIO->next_seq;
    my @organism = $seqobj->species->classification;
    my $strain   = $organism[$tax];
    my %tmplist;
    $tmplist{$strain} = $i;
    $pm->finish( 0, \%tmplist );
  }
  return %genomelist;
}

###################################################################################
sub GenerateKmers() {
###################################################################################
  my ( $in, $outfile, $informat, $k, $taxid ) = @_;
  my %oligo;
  my @files = @$in;
  foreach my $infile (@files) {
    my $seqIO = Bio::SeqIO->new( -file => "$infile", -format => "$informat" );
    while ( my $seqobj = $seqIO->next_seq ) {
      my $seq = $seqobj->seq;
      for ( my $i = 0; $i <= ( length($seq) - $k - 1 ); $i++ ) {
        my $oligo = substr( $seq, $i, $k );
        $oligo{$oligo}++ if $oligo !~ /[^ATCGatcg]/;
      }
    }
  }
  my $j = 1;
  open( FA, ">$outfile" ) || die "# can not write $outfile in Generate k-mers\n";
  foreach my $oligo ( keys %oligo ) {
    print FA ">$taxid\_$j\n$oligo\n";
    $j++;
  }
  close FA;
}

###################################################################################
sub GenerateKmerDB() {
###################################################################################
  my ( $infafile, $kmerdbfile, $meryl, $segment, $threads, $k ) = @_;
  system(
    "$meryl -B -v -C -m $k -segments $segment -threads $threads -s $infafile -o $kmerdbfile 1>/dev/null 2>/dev/null"
  );
}

###################################################################################
sub GenerateGSMs() {
###################################################################################
  my ( $infile, $outfile, $informat, $k, $taxid ) = @_;
  my %oligo;
  my $seqIO = Bio::SeqIO->new( -file => "$infile", -format => "$informat" );
  while ( my $seqobj = $seqIO->next_seq ) {
    my $seq = $seqobj->seq;
    for ( my $i = 0; $i <= ( length($seq) - $k - 1 ); $i++ ) {
      my $oligo = substr( $seq, $i, $k );
      $oligo{$oligo}++ if ( $oligo !~ m/[^ATCGatcg]/ );
    }
  }
  my $j = 1;
  open( FA, ">$outfile" ) || die "# can not write $outfile in Generate k-mers\n";
  foreach my $oligo ( keys %oligo ) {
    print FA ">$taxid\-$j\n$oligo\n";
    $j++;
  }
  close FA;
}

###################################################################################
sub MapAndFilter() {
###################################################################################
  my ( $kmerdb, $infafile, $outmapfile, $outfafile, $k ) = @_;
  system("$mapMers -m $k -mers $kmerdb -seq $infafile > $outmapfile  2>/dev/null");

  #system("$mapMers -m $k -mers $kmerdb -seq $infafile > $outmapfile");
  ##read map file
  my %map;
  open( MAP, "$outmapfile" ) || die "#1\n";
  while (<MAP>) {
    chomp;
    my @items = split( "\t", $_ );
    $map{ $items[0] } = 1;
  }
  close MAP;
  ##filter mapped probe mers
  $/ = ">";
  open( FASTA, "$infafile" )   || die "can't open file\n";
  open( OUT,   ">$outfafile" ) || die "can't open file\n";
  my $junk = (<FASTA>);
  while ( my $frecord = <FASTA> ) {
    chomp $frecord;
    my ( $fdef, @seqLines ) = split /\n/, $frecord;
    my $seq = join( "", @seqLines );
    print OUT ">$fdef\n$seq\n" if !$map{$fdef};
  }
  close OUT;
  close FASTA;
  $/ = "\n";
}

###################################################################################
sub CheckProbeSpecificity() {
###################################################################################
  my ( $probefile, $blastfile, $seqfile, $probeoutfile ) = @_;
  my %probes2keep;
  my %nontarget;
  my %probes = &GetProbeSeq($probefile);
  my %blastn = &ReadBLASTN( $blastfile, $tax );
  my ( $nontarget_ids, $target_ids ) = &FetchSeqIDFromBLASTN( \%blastn, $tax );
  my @nontarget_ids = @$nontarget_ids;
  my @target_ids    = @$target_ids;

  #my %gseq          = &GetGenomeSeq( $seqfile, \@nontarget_ids );
  my %numcds;
  my %allgeneinfo;
  foreach my $gid (@target_ids) {
    my $gbkfile  = "$gbkdir/${gid}.gbk";
    my %geneinfo = &GetGeneInfo($gbkfile);
    %allgeneinfo = ( %allgeneinfo, %geneinfo );
    my $numcds = `grep -c "^     CDS" gbk/$gid.gbk`;
    chomp $numcds;
    $numcds{$gid} = $numcds;
  }

  open( PROBEOUT, ">>$probeoutfile" ) || die "# can not write $probeoutfile\n";
  foreach my $probeID ( keys %probes ) {
    my $probe = $probes{$probeID};
    $probeID =~ /(\d+)\_(\d+)\-\d+/;
    my ( $i, $j ) = ( $1, $2 );
    my $genomeID = $i if $tax == 1;
    $genomeID = "$i\_$j" if $tax == 0;
    my @positions;
    my @start;
    my @end;

    if ( $blastn{$probeID} ) {
      foreach my $hitID ( keys %{ $blastn{$probeID} } ) {
        if ( $hitID !~ /^$genomeID(-|_)/ ) {

          my $gidentity =
            $blastn{$probeID}{$hitID}{aln} *
            $blastn{$probeID}{$hitID}{identity} /
            $probe_length;

          if ( !$nontarget{$probeID} ) {
            $nontarget{$probeID}{gidentity} = $gidentity;

          }
          else {
            if ( $gidentity > $nontarget{$probeID}{gidentity} ) {
              $nontarget{$probeID}{gidentity} = $gidentity;

            }
          }
        }
        else {
          if ( $blastn{$probeID}{$hitID}{aln} == $probe_length
            && $blastn{$probeID}{$hitID}{identity} == 100 )
          {
            $hitID =~ /^($genomeID)(.*?)\-(.*?)$/;
            my $seqid  = $3;
            my $gid    = "$1" . "$2";
            my %gene   = %{ $allgeneinfo{$seqid} } if $allgeneinfo{$seqid};
            my $hstart = $blastn{$probeID}{$hitID}{hstart};
            my $hend   = $blastn{$probeID}{$hitID}{hend};
            push( @start, "$hitID\_$hstart" );
            push( @end,   "$hitID\_$hend" );

            if ( $hstart < $hend ) {
              my $position = &IdentifyProbePos( $hstart, $hend, \%gene );
              if ( $numcds{$gid} > 0 ) {
                push( @positions, $position );
              }
              else {
                push( @positions, "unclear" );
              }
            }
            else {
              my $position = &IdentifyProbePos( $hend, $hstart, \%gene );
              if ( $numcds{$gid} > 0 ) {
                push( @positions, $position );
              }
              else {
                push( @positions, "unclear" );
              }
            }
          }
        }
      }
      my %temp;
      @positions = grep { !$temp{$_}++ } @positions;
      if ( $nontarget{$probeID} ) {
        my $gidentity = $nontarget{$probeID}{gidentity};

        #my $stretch   = $nontarget{$probeID}{stretch};
        my $length = length($probe);
        my $start  = join( ";", @start );
        my $end    = join( ";", @end );

        if ( $gidentity <= $identity_cutoff ) {
          print PROBEOUT "$probeID\t$length\t$start\t$end\t$gidentity\t@positions\t$probe\n";
        }
      }
      else {
        my $probescore = 1;
        my $length     = length($probe);
        my $start      = join( ";", @start );
        my $end        = join( ";", @end );

        if ( $length == $probe_length ) {
          print PROBEOUT "$probeID\t$length\t$start\t$end\tna\t@positions\t$probe\n";
        }
      }
    }
    else {
      print "ERROR:$probeID has no hit with database!\n";
    }
  }
  close PROBEOUT;
}

###################################################################################
sub GetProbeSeq() {
###################################################################################
  my ($file) = @_;
  my %probes;
  $/ = ">";
  open( FASTA, "$file" ) || die "# can't open genome sequence file $file\n";
  my $junk = <FASTA>;
  while ( my $line = <FASTA> ) {
    chomp $line;
    my ( $idline, @seqlines ) = split /\n/, $line;
    my $seq = join '', @seqlines;
    my @token = split( " ", $idline );
    $probes{ $token[0] } = $seq;
  }
  close FASTA;
  $/ = "\n";
  return %probes;
}

###################################################################################
sub GetGenomeSeq() {
###################################################################################
  my ( $file, $array ) = @_;
  my ( @array, %array );
  if ($array) {
    @array = @$array;
    %array;
    foreach my $id (@array) {
      $array{$id} = 1;
    }
  }
  my %genomeseq;
  $/ = ">";
  open( FASTA, "$file" ) || die "# can't open genome sequence file $file\n";
  my $junk = <FASTA>;
  while ( my $line = <FASTA> ) {
    chomp $line;
    my ( $idline, @seqlines ) = split /\n/, $line;
    my $seq = join '', @seqlines;
    my @token = split( " ", $idline );
    if ( $#array >= 0 ) {
      if ( $array{ $token[0] } ) {
        $genomeseq{ $token[0] } = $seq;
      }
    }
    else {
      $genomeseq{ $token[0] } = $seq;
    }
  }
  close FASTA;
  $/ = "\n";
  return %genomeseq;
}

###################################################################################
sub ReadBLASTN() {
###################################################################################
  my ( $file, $tax ) = @_;
  my %blastn;
  open( BLAST, "$file" ) || die "# can not open blast file $file\n";
  while (<BLAST>) {
    chomp;
    my @items = split( "\t", $_ );
    $items[0] =~ /(\d+)\_(\d+)\-\d+/;
    my ( $i, $j ) = ( $1, $2 );
    my $gid = $i if $tax == 1;
    $gid = "$i\_$j" if $tax == 0;
    if (
        !$blastn{ $items[0] }{ $items[1] }
      && $items[3] >= $stretch_cutoff
      && ( ( $items[1] =~ /^$gid(-|_)/ )
        || ( $items[1] !~ /^$gid(-|_)/ && ( $items[2] < 100 || $items[3] < $probe_length ) ) )
      )
    {
      $blastn{ $items[0] }{ $items[1] }{qstart}   = $items[6];
      $blastn{ $items[0] }{ $items[1] }{qend}     = $items[7];
      $blastn{ $items[0] }{ $items[1] }{hstart}   = $items[8];
      $blastn{ $items[0] }{ $items[1] }{hend}     = $items[9];
      $blastn{ $items[0] }{ $items[1] }{aln}      = $items[3];
      $blastn{ $items[0] }{ $items[1] }{identity} = $items[2];
    }
  }
  close BLAST;
  return %blastn;
}

###################################################################################
sub FetchSeqIDFromBLASTN() {
###################################################################################
  my ( $blastn, $tax ) = @_;
  my %blastn = %$blastn;
  my ( @ids1, @ids2 );    ##non target genomeIDs, and target genomeIDs
  foreach my $probeID ( keys %blastn ) {
    $probeID =~ /(\d+)\_(\d+)\-\d+/;
    my $i        = $1;
    my $j        = $2;
    my $genomeID = $i if $tax == 1;
    $genomeID = "$i\_$j" if $tax == 0;

    #push( @ids2, "$i\_$j" );
    foreach my $hit ( keys %{ $blastn{$probeID} } ) {
      if ( $hit !~ /^$genomeID(_|-)/ ) {
        push( @ids1, $hit );
      }
      else {
        if ( $hit =~ /^(\d+)\_(\d+)\-/ ) {
          push( @ids2, "$1\_$2" );
        }
      }
    }
  }
  my %temp;
  @ids1 = grep { !$temp{$_}++ } @ids1;
  %temp = "";
  @ids2 = grep { !$temp{$_}++ } @ids2;
  return ( \@ids1, \@ids2 );
}

###################################################################################
sub GetTargetRegion() {
###################################################################################
  my ( $hseq, $qstart, $qend, $hstart, $hend, $probelength ) = @_;
  my $seq;
  if ( $hend > $hstart ) {
    my $startpos = $hstart - $qstart;
    $seq = substr( $hseq, $startpos, $probelength );
  }
  else {
    my $startpos = $hend - ( $probelength - $qend + 1 );
    $seq = substr( $hseq, $startpos, $probelength );
    $seq = reverse($seq);
    $seq =~ tr/ATCGatcg/TAGCtagc/;
  }
  return $seq;
}

###################################################################################
sub CalStretch() {
###################################################################################
  my ( $seq1, $seq2 ) = @_;
  my @n1 = split( "", $seq1 );
  my @n2 = split( "", $seq2 );
  my @stretch;
  my $stretch;
  for ( my $length = 50; $length >= 1; $length-- ) {
    for ( my $start = 0; $start <= 50 - $length; $start++ ) {
      my $substr1 = substr( $seq1, $start, $length );
      my $substr2 = substr( $seq2, $start, $length );
      if ( $substr1 eq $substr2 ) {
        push( @stretch, $length );
      }
    }
  }
  @stretch = sort { $a <=> $b } @stretch;
  $stretch = $stretch[$#stretch];
  return $stretch;
}

###################################################################################
sub GetGeneInfo() {
###################################################################################
  my ($gbkfile) = @_;
  my %gene;
  my %annotation;
  my $seqio = Bio::SeqIO->new( -file => "$gbkfile", -format => "genbank" );
  while ( my $seqobj = $seqio->next_seq ) {
    my $seqid    = $seqobj->id;
    my @features = $seqobj->get_SeqFeatures();
    for ( my $i = 0; $i < @features; $i++ ) {
      my $feat = $features[$i];
      if ( $feat->primary_tag eq "CDS" ) {
        my %feat       = %$feat;
        my $cds_start  = $feat{_location}{_start};
        my $cds_end    = $feat{_location}{_end};
        my $cds_strand = $feat{_location}{_strand};
        my @cds_locus_tag;
        if ( $feat->has_tag("locus_tag") ) {
          @cds_locus_tag = $feat->get_tag_values("locus_tag");
        }
        else { @cds_locus_tag = ("$cds_start\_$cds_end"); }
        my @cds_product = $feat->get_tag_values("product") if $feat->has_tag("product");
        $annotation{ $cds_locus_tag[0] } = $cds_product[0];
      }
      if ( $feat->primary_tag eq "gene" ) {
        my %feat        = %$feat;
        my $gene_start  = $feat{_location}{_start};
        my $gene_end    = $feat{_location}{_end};
        my $gene_strand = $feat{_location}{_strand};
        my @gene_locus_tag;
        if ( $feat->has_tag("locus_tag") ) {
          @gene_locus_tag = $feat->get_tag_values("locus_tag");
        }
        else { @gene_locus_tag = ("$gene_start\_$gene_end"); }
        $gene{$seqid}{ $gene_locus_tag[0] }{start}  = $gene_start;
        $gene{$seqid}{ $gene_locus_tag[0] }{end}    = $gene_end;
        $gene{$seqid}{ $gene_locus_tag[0] }{strand} = $gene_strand;
        $gene{$seqid}{ $gene_locus_tag[0] }{type}   = "gene";
      }
    }
  }
  foreach my $seqid ( keys %gene ) {
    foreach my $locus_tag ( keys %{ $gene{$seqid} } ) {
      if ( $annotation{$locus_tag} ) {
        $gene{$seqid}{$locus_tag}{product} = $annotation{$locus_tag};
      }
    }
  }
  return %gene;
}

###################################################################################
sub IdentifyProbePos() {
###################################################################################
  my ( $start, $end, $geneinfo ) = @_;
  my %geneinfo = %$geneinfo;
  my %position;
  my $position;
  foreach my $locus_tag ( keys %geneinfo ) {
    if ( $start >= $geneinfo{$locus_tag}{start} && $end <= $geneinfo{$locus_tag}{end} ) {
      $position{"gene"} = 1;
    }
    if (
      (
           $start >= $geneinfo{$locus_tag}{start}
        && $start <= $geneinfo{$locus_tag}{end}
        && $end > $geneinfo{$locus_tag}{end}
      )
      || ( $start <= $geneinfo{$locus_tag}{start}
        && $end >= $geneinfo{$locus_tag}{start}
        && $end <= $geneinfo{$locus_tag}{end} )
      )
    {
      $position{"pgene"} = 1;
    }
  }
  my @position = keys %position;
  if ( $#position == -1 ) {
    $position = "intergenic";
  }
  else {
    if ( $#position == 1 ) {
      $position = "gene";
    }
    else {
      $position = $position[0];
    }
  }
  return $position;
}

###################################################################################
sub printhelp() {
###################################################################################
  print STDERR <<ENDHELP;
Parameters of $0:
perl GSMer.pl -m ("splitgbk"|"readgbk"|"makeblastdb"|"makekmerdb"|"getgsm"|"mapgsm"|"blastgsm"|"checkspecificity") -i <input merged genbank file, required for "splitgbk"> -l <list of interested microbial strains, optional> -f1 <input alien fasta file, e.g. human genome>
Tutorials: 
    -m: different modules, need to run step by step, available options are:
        step 1a. "splitgbk"--split huge merged genbank file into strain/species level gbk files, a strain.list file will be generated
                             e.g. perl GSMer.pl -m splitgbk -i all.gbk
        step 1b. "readgbk"--read gbk file from gbk directory and generate a strain.list file. This is an alternative option when strain/species level gbk files are already available.
                            e.g. perl GSMer.pl -m readgbk
        step 2. "makeblastdb"--make blast database from all microbial genomes and alien genomes if available
                            e.g. perl GSMer.pl -m makeblastdb -f1 human.fasta
        step 3. "makekmerdb"--generate k-mer database for k-mer mapping, only k-mers showing up in >=2 genomes were used; if alien genome is provided, all k-mers from alien genome will be used.                    
                            e.g. perl GSMer.pl -m makekmerdb -f1 human.fasta
        step 4. "getgsm"--generate all candidate GSMs from microbial genomes, e.g. to get all GSMs for e.coli strains: 
                          e.g. perl GSMer.pl -m getgsm -l ecoli.list  
        step 5. "mapgsm"--map GSMs to k-mer database forcontinuous stretch filtering
                          e.g. perl GSMer.pl -m mapgsm -l ecoli.list
        step 6. "blastgsm"--blast search GSMs against blast database for identity filtering
                          e.g. perl GSMer.pl -m blastgsm -l ecoli.list   
        step 7. "checkspecificity"--filtering blast results, *.out files are generated for each strain/species
                          e.g. perl GSMer.pl -m checkspecificity -l ecoli.list
ENDHELP
  exit;
}

####################################################################
sub print_start(){
####################################################################
print STDOUT <<START;
------------------------------------------------------------------------------------------------------------
					Welcome to GSMer!
------------------------------------------------------------------------------------------------------------
START
}


####################################################################
sub print_end(){
####################################################################
print STDOUT <<END;
------------------------------------------------------------------------------------------------------------
					Thanks for using GSMer!
------------------------------------------------------------------------------------------------------------      
END
}
