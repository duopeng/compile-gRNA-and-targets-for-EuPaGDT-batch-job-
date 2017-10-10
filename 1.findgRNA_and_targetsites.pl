#! /usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage; 
my $genome='';
my $destdir='';
my $message_text  = "please specify option parameters. -i < job folder from EuPaGDT> -g <genome, specified in the EuPaGDT job>";
  my $exit_status   = 2;          ## The exit status to use
  my $verbose_level = 0;          ## The verbose level to use
  my $filehandle    = \*STDERR;   ## The filehandle to write to
GetOptions ('i=s' => \$destdir, 'g=s' => \$genome);
pod2usage( { -message => $message_text ,
               -exitval => $exit_status  ,  
               -verbose => $verbose_level,  
               -output  => $filehandle } ) if ($destdir eq '' or $genome eq '');
$destdir=$destdir.'/genes';

my @fields=split(/\//,$destdir);
my $jobname=$fields[$#fields-1];

#####################################################
#go through the EuPaGDT result
#####################################################
my %TcTSgRNAsummaryoff;
my %TcTSgRNAsummaryon;
my %TcTSgRNA_newTcTScollection_ontarget;
my $currentgRNA='';
my %gRNAseq; # $gRNAseq{$gRNA}=$seq
opendir DIR, "$destdir" or die;
while(my $filename=readdir(DIR))
{  
	next if $filename eq "." or $filename eq "..";
 #######on target
	open ON, "$destdir/${filename}/${filename}.gRNA.fasta.${genome}.blastout.parsed.on.html" or next;
	<ON>;<ON>;<ON>; # header
	while (my $line=<ON>)
	{
		chomp $line;		
		if ($line=~/colspan\=1/) # gRNA name&seq line
		{
			$line=~s/<tr><td colspan=1>//g;
			$line=~s/<\/td><td colspan=7><\/td><\/tr>//g;
			$line=~s/<\/td><\/tr>//;
			$line=~s/ colspan=8//g;
			my ($gRNAname,$gRNAseq)=split (/<\/td><td>/, $line);
			$currentgRNA=$gRNAname;
			$gRNAseq{$gRNAname}=$gRNAseq;
		}
		if ($line=~/colspan\=2/) # gRNA info line
		{
			$line=~s/<tr><td colspan=2>//g;
			$line=~s/<label style="color:#FF0000;font-weight:bold">//g;
			$line=~s/<label style="color:#FF9900;font-weight:bold">//g;
		    $line=~s/<\/td><\/tr>//g;
			$line=~s/<\/label>//g;
			$line=~s/ colspan=8//g;
			my ($empty,$gRNA_match_start,$gRNA_match_end,$alignment,$ontargettext,$alignment_identity,$genome_annotation,$match_chromosome,$chromosome_match_start,$chromosome_match_end)=split (/<\/td><td>/, $line);
			if (defined $gRNA_match_end)
			{
				$gRNA_match_end=~s/<\/td>//;
			if (defined $chromosome_match_end)
		  {
			if ($gRNA_match_end==27 and $alignment_identity>=90)  # is on target
			{
				if ((exists $TcTSgRNAsummaryon{$currentgRNA}) ) 
					{
						$TcTSgRNAsummaryon{$currentgRNA}+=1;
					}
					else {
						$TcTSgRNAsummaryon{$currentgRNA}=1;
					}
			}			
		 }
		 }
		}	
	}
	close ON;
#######off target
	open OFF, "$destdir/${filename}/${filename}.gRNA.fasta.${genome}.blastout.parsed.off.html" or next;
	<OFF>;<OFF>;<OFF>; # header
	while (my $line=<OFF>)
	{
		chomp $line;	
		if ($line=~/colspan\=1/) # gRNA name&seq line
		{
			$line=~s/<tr><td colspan=1>//g;
			$line=~s/<\/td><td colspan=7><\/td><\/tr>//g;
			$line=~s/<\/td><\/tr>//;
			$line=~s/ colspan=8//;
			
			my ($gRNAname,$gRNAseq)=split (/<\/td><td>/, $line);
			$currentgRNA=$gRNAname;
			$TcTSgRNAsummaryoff{$currentgRNA}=0;
		}
		if ($line=~/colspan\=2/) # gRNA info line
		{			
			if ($line=~/no off-target hits/) #no offtarget
			{				
				#print "no off target\n"
			}
			else {	#has off target			
				$line=~s/<tr><td colspan=2>//g;
				$line=~s/<label style="color:#FF0000;font-weight:bold">//g;
				$line=~s/<label style="color:#FF9900;font-weight:bold">//g;
				$line=~s/<\/td><\/tr>//g;
				$line=~s/<\/label>//g;	
				$line=~s/ colspan=8//g;
				my ($empty,$gRNA_match_start,$gRNA_match_end,$alignment,$ontargettext,$alignment_identity,$genome_annotation,$match_chromosome,$chromosome_match_start,$chromosome_match_end)=split (/<\/td><td>/, $line);
				if (defined $gRNA_match_end)
				{
	               $gRNA_match_end=~s/<\/td>//;
				if (defined $chromosome_match_end)
				{
					if ($gRNA_match_end==27 and $alignment_identity>=85 and !($line=~/trans-sialidase/))  # is off target
					{			
						my ($geneid,$anno)=split(':',$genome_annotation);			
						if (exists $TcTSgRNAsummaryoff{$currentgRNA})
						{
						$TcTSgRNAsummaryoff{$currentgRNA}+=1;
						}
						else {
							$TcTSgRNAsummaryoff{$currentgRNA}=1;
						}			
				}			
				}
				}
			}
		}		
	}
	close OFF; 
}

##############################
####### summarize gRNA targets
##############################
## on targets
open ONOUT, ">${jobname}.on_target_number_list.genome_off_included.txt";
print ONOUT "gRNA\tseq\ton-target\toff-target\n";
foreach my $gRNA (sort keys %TcTSgRNAsummaryon)  # deal with zeros
{	
	$TcTSgRNAsummaryon{$gRNA}=0 if $TcTSgRNAsummaryon{$gRNA} eq "";;
	if ($gRNA eq "")
	{
		delete $TcTSgRNAsummaryon{$gRNA};
		#print "deleted $gRNA\n"
	}
}
foreach my $gRNA (sort keys %TcTSgRNAsummaryoff) # deal with zeros
{
	
	$TcTSgRNAsummaryoff{$gRNA}=0 if $TcTSgRNAsummaryoff{$gRNA} eq "";;
	if ($gRNA eq "")
	{
		delete $TcTSgRNAsummaryoff{$gRNA};
		#print "deleted $gRNA\n"
	}
}
foreach my $gRNA (sort {$TcTSgRNAsummaryon{$b} <=> $TcTSgRNAsummaryon{$a} } keys %TcTSgRNAsummaryon) #print result
{
	my $num=$TcTSgRNAsummaryon{$gRNA};
	my $offnum=0;
	$offnum=$TcTSgRNAsummaryoff{$gRNA} if (exists $TcTSgRNAsummaryoff{$gRNA});	
	my $gRNAseq=$gRNAseq{$gRNA};	
	my $newTcTS_collection_on=$TcTSgRNA_newTcTScollection_ontarget{$gRNA};	
	print ONOUT "$gRNA\t$gRNAseq\t$num\t$offnum\n";
}
close ONOUT;

