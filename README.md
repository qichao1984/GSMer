GSMer 
=====
GSMer devotes to identify genome-specific markers (GSMs) from currently sequenced microbial genomes using a k-mer based approach. The GSMs could be used to identify microbial strains/species in metagenomes, especially in human microbiome where many reference genomes are available. Two different levels of GSMs, including strain-specific and species-specifc GSMs are currently supported. 

Dependencies:  
-----------------------------------------------------------------------
<dl>
<dt>Third party programs:</dt>
</dl>
* NCBI BLAST (megablast+formatdb)  
ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.26/ncbi.tar.gz  
* KMER (meryl+mapMers)  
http://kmer.sourceforge.net

<dl>
<dt>Perl libraries:</dt>
</dl>
* Bio::SeqIO  
* Getopt::Long  
* Parallel:ForkManager  
* String::Random  

Tutorial:
-----------------------------------------------------------------------
