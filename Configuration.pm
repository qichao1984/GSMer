#
#===============================================================================
#
#         FILE: Configuration.pm
#
#  DESCRIPTION: Configuration file for GSM identification
#
#        FILES: Configuration.pm
#         BUGS: ---
#        NOTES: set program path and parameters here for GSM identification
#       AUTHOR: Qichao Tu (qichao@ou.edu, philloid@gmail.com)
#  INSTITUTION: Institute for Environmental Genomics, University of Oklahoma
#      VERSION: 1.0
#      CREATED: 11/16/2013 5:11:38 PM
#     REVISION: ---
#===============================================================================
package Configuration;
{
####
  # SYSTEM PARAMETERS
  #
  # ENTER HERE THE PATH OF THE PROGRAMS
  # Make sure you do not accidently delete the semicolon at the end of the lines
  ## level of specificity, 2 for species level, 1 for strain level
  $tax=1;
#####
  ##blastall exe file, blast search for similar non-targets
  $blastall = "./exe/megablast";
  ##formatdb exe file
  $formatdb = "./exe/formatdb";
  ##meryl exe file, k-mer database generation, extraction, and modification
  $meryl = "/opt/apps/kmer/bin/meryl";
  ##mapMers exe file, map k-mers to probe mers
  $mapMers = "/opt/apps/kmer/bin/mapMers";
  #####
  # PROJECT PARAMETERS
  #####

  ##dir for k-mer databases, gbk files, output dir
  #k-mer database directory
  $merdbdir = "merdb_strain";
  #gbk file directory
  $gbkdir   = "gbk";
  #marker output directory
  $probedir   = "probe_species_gut";
  #blast database directory
  $blastdbdir   = "blastdb";
  ## ROOT DIRECTORY FOR ANALYSIS 
  $root_directory = ".";
  ## multi-core support
  $threads = "35";
  $segment=1;#segment number for meryl to speed up
  ##probe length, stretch, identity etc
  #k-mer size for stretch filtering
  $k                     = 18; 
  #marker length
  $probe_length          = 50;
  #continuous stretch cutoff
  $stretch_cutoff        = 20;
  #sequence identity cutoff
  $identity_cutoff       = 85;
}

1;
