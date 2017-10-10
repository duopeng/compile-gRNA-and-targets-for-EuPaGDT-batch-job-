#! /usr/bin/perl
# used to collapse same gRNA in the list
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Perl;
use Getopt::Long;
use Pod::Usage; 
my $infile="";
my $fasta='';
my $message_text  = "please specify option parameters. -i < output file from 1.findgRNA_and_targetsites.pl script > -f <fasta file of all genes used to run the EuPaGDT job, for confirming on-targets>";
  my $exit_status   = 2;          ## The exit status to use
  my $verbose_level = 0;          ## The verbose level to use
  my $filehandle    = \*STDERR;   ## The filehandle to write to
GetOptions ('i=s' => \$infile, 'f=s' =>\$fasta);
pod2usage( { -message => $message_text ,
               -exitval => $exit_status  ,  
               -verbose => $verbose_level,  
               -output  => $filehandle } ) if ($infile eq '' or $fasta eq '');
			   
#######read in all TcTS
my %allTcTS;
my $seqio=Bio::SeqIO->new(-format => 'Fasta', -file=>"$fasta");
while(my $seqobj=$seqio->next_seq)
{
	my $id=$seqobj->id;
	my $seq=$seqobj->seq;
	$allTcTS{$id}=$seq;
	#print "$seq\n"
}

my %gRNA;			   
my %sort_gRNA;
open IN, "$infile";
my $header=<IN>; #header
while (my $line=<IN>)
{
	chomp $line;
	my ($gRNA,$seq,$on_target,$off_target,$on_target_new)=split("\t",$line);
	if (exists $gRNA{$seq})
	{
		if ($gRNA{$seq}{'$on_target'}<$on_target)
		{
			$gRNA{$seq}{'line'}=$line;
			$gRNA{$seq}{'$on_target'}=$on_target;
			$sort_gRNA{$seq}=$on_target;
		}
	}
	else
	{
		$gRNA{$seq}{'line'}=$line;
		$gRNA{$seq}{'$on_target'}=$on_target;
		$sort_gRNA{$seq}=$on_target;
	}
}

close IN;
open OUT, ">${infile}.collapsed";

chomp $header;
print OUT $header."on-target(in fasta file)\n";
foreach my $key (sort {$sort_gRNA{$b} <=> $sort_gRNA{$a} } keys %sort_gRNA)
{
			my $newTcTS_ontarget=0;
            my $gRNAseq=$key;
			$gRNAseq=~s/\|//;
			my $revcom_seqobj = revcom( $gRNAseq );
			my $gRNAseq_revcom=$revcom_seqobj->seq();
			##find on target in new defined TcTS collection
			foreach my $keys (sort keys %allTcTS)
			{
				if ($allTcTS{$keys}=~/$gRNAseq/i or $allTcTS{$keys}=~/$gRNAseq_revcom/i )
				{	
					  $newTcTS_ontarget+=1;
				}	
			}
	print OUT $gRNA{$key}{'line'}."\t$newTcTS_ontarget\n";
}