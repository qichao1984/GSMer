GSMer 
=====
GSMer devotes to identify genome-specific markers (GSMs) from currently sequenced microbial genomes using a k-mer based approach. Explored GSMs could be used to identify microbial strains/species in metagenomes, especially in human microbiome where many reference genomes are available. Two different levels of GSMs, including strain-specific and species-specifc GSMs are currently supported. 

![](https://github.com/qichao1984/GSMer/blob/master/GSMer_pipeline.png =300x250 "GSMer pipeline")

####Citation:
#####Qichao Tu, Zhili He and Jizhong Zhou. “Strain/Species identification in metagenomes using genome-specific markers.” Nucleic Acids Research. (accepted)

Identified GSMs
=====
Species-specific GSMs: 2,005 species (4,933 strains). 
* All species-specific GSMs: 
http://ieg.ou.edu/GSMer/allgsm_species.zip
* Randomly selected upto 100 GSMs/strain from different regions:
http://ieg.ou.edu/GSMer/gsm_species.select100.species.zip

Strain-specific GSMs: 4,088 strains
* All strain-specific GSMs:
http://ieg.ou.edu/GSMer/allgsm_strain.zip
* Randomly selected upto 50 GSMs/strain from different regions: 
http://ieg.ou.edu/GSMer/gsm_strain.select50.strains.zip

Full list of included microbial strains:
* http://ieg.ou.edu/GSMer/strain.xlsx

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
* Bio::SeqIO (bioperl)  
* Getopt::Long  
* Parallel:ForkManager  
* String::Random  

Tutorial:
-----------------------------------------------------------------------
This tutorial shows how to identify GSMs for E.coli O157 with E.coli K12 genome as alien.  
A total of seven steps are required for GSM identification, and need to run one by one. For details of all available options, please run *perl GSMer.pl -m help*.  
<dl>
<dt>Testing files:</dt>
</dl>
* O157.gbk: E.coli O157 genomes in genbank format, four differet strains were included (i.e. O157:H7 EC4115, O157:H7 EDL933, O157:H7 TW14359, and O157:H7 Sakai). 
* k12.fa: E.coli K-12 substr. W3110 genome in fatsa format. This genome will be used alien genome, resulting in O157-specific GSMs, which would not be found in the K-12 genome.

<dl>
<dt>Steps to run:</dt>
</dl>
0. Check the Configuration.pm file, and set tax level at 1, which represent strain level. Make sure all other program path and GSM criteria are correct. 
1. `perl GSMer.pl -m splitgbk -i O157.gbk`  
   This step split the O157.gbk file into four gbk files representing the four O157 strains. Four gbk files will be generated in a gbk directory. A strain.list file will also be generated in the working directory. 
2. `GSMer.pl -m makeblastdb -f1 k12.fa`  
   This step create a blast database file from all the four gbk files in the gbk directory, as well as the K-12 genome. 
3. `perl GSMer.pl -m makekmerdb -f1 k12.fa`  
   This step create a k-mer database for all O157 genomes and k-12 genome. K-mers that show up in >=2 O157 genomes and all k-mers in K-12 genome are extracted for k-mer database construction.  
4. `perl GSMer.pl -m getgsm`  
   This step generate all candidate GSMs for O157 strains.  
5. `perl GSMer.pl -m mapgsm`  
   This step maps the above candidate GSMs to the k-mer database for continuous stretch filtering.  
6. `perl GSMer.pl -m blastgsm`  
   This step performs blast searching unmapped GSMs against the blast database.  
7. `perl GSMer.pl -m checkspecificity`  
   This step filters GSMs based the blast output for continuous stretch and identity with non-target genomes. Four *.out files containing detailed information of O157-specific GSMs will be generated in the GSM output directory specified in the Configuration.pm file.
