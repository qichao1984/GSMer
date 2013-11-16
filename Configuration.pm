#
#===============================================================================
#
#         FILE: Configuration.pm
#
#  DESCRIPTION:
#
#        FILES: Configuration.pm
#         BUGS: ---
#        NOTES: set program path and parameters here for GSM identification
#       AUTHOR: Qichao Tu
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
  ## level of specificity, 1 for species level, 0 for strain level
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
  $merdbdir = "merdb_strain";
  $gbkdir   = "gbk";
  $probedir   = "probe_species_gut";#marker output dir
  $blastdbdir   = "blastdb";
  ## ROOT DIRECTORY FOR ANALYSIS (absolute pathname required)
  $root_directory = ".";
  ## multi-core
  $threads = "35";
  ##probe length, stretch, identity etc
  $k                     = 18;
  $probe_length          = 50;
  $stretch_cutoff        = 20;
  $identity_cutoff       = 90;
}

1;
