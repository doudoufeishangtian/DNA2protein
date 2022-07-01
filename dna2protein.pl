#!/usr/bin/perl
=head1 Command-line Option
  
  perl cds2aa.pl  <infile.fa | STDIN>
  --translate <str>   set the translate table, default=standard
  --check             check the quality of gene model
  --verbose           output verbose information to screen  
  --help              output help information to screen  

=head1 Usage Exmples

  perl ./cds2aa.pl test-data/BMF000017.cds.fa 
  perl ./cds2aa.pl test-data/BMF000017.exon.fa 
  perl ./cds2aa.pl -check test-data/test_3chrs.fa.fgenesh.100.cds 

=cut    #注释结束

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Data::Dumper;


my ($Check,$Trans_table,$Verbose,$Help);
GetOptions(
    "translate:s"=>\$Trans_table,
    "check"=>\$Check,
    "verbose"=>\$Verbose,
    "help"=>\$Help
) || die "Please use --help option to get help\n";
$Trans_table ||= "standard";
die `pod2text $0` if ($Help);


my ($Check_start,$Check_stop,$Check_mid,$Check_triple) = (0,0,0,0);

my %CODE = (
            "standard" =>
                {   
                'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
                'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
                'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
                'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
                'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
                'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
                'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
                'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
                'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
                'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
                'ATG' => 'M',                                                                         # Methionine
                'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
                'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
                'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
                'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
                'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
                'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
                'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
                'TGG' => 'W',                                                                         # Tryptophan
                'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
                'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
                }
            ## more translate table could be added here in future
            ## more translate table could be added here in future
            ## more translate table could be added here in future
    );


print "Id\tstart\tstop\tmiddle\ttriple\n" if($Check);

$/=">"; <>; $/="\n";
while (<>) {
    my $head = $_;
    chomp $head;
    my $key = $1 if($head =~ /^(\S+)/);
    my $phase = ($head =~ /\s+phase[:\s]+([012])\s+/i) ? $1 : 0 ;
    $/=">";
    my $seq = <>;
    chomp $seq;
    $/="\n";
    
    if ($Check) {
        my ($start,$end,$mid,$triple) = check_CDS($seq);
        print "$key\t$start\t$end\t$mid\t$triple\n" if(!$start || !$end || !$mid || !$triple);
        $Check_start++   if($start == 0);
        $Check_stop++   if($end == 0);
        $Check_mid++   if($mid == 0);
        $Check_triple++  if($triple == 0);
    }else{
        my $prot = cds2aa($seq,$phase,$CODE{$Trans_table});
        Display_seq(\$prot);
        print ">$head [translate_table: $Trans_table]\n".$prot;
    }

}
close IN;

if ($Check && $Verbose){
    print STDERR "wrong_start\twrong_stop\twrong_middle\twrong_triple\n";
    print STDERR "$Check_start\t$Check_stop\t$Check_mid\t$Check_triple\n";
}


####################################################
################### Sub Routines ###################
####################################################


#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#       disp_seq(\$string);
#############################################
sub Display_seq{
    my $seq_p=shift;
    my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
    my $disp;

    $$seq_p =~ s/\s//g;
    for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
        $disp .= substr($$seq_p,$i,$num_line)."\n";
    }
    $$seq_p = ($disp) ?  $disp : "\n";
}
#############################################


## translate CDS to pep
####################################################
sub cds2aa {
    my $seq = shift;
    my $phase = shift || 0; 
    my $translate_p = shift;
    
    $seq =~ s/\s//g;
    $seq = uc($seq);
    
    my $len = length($seq);
    
    my $prot;
    for (my $i=$phase; $i<$len; $i+=3) {
        my $codon = substr($seq,$i,3);
        last if(length($codon) < 3);
        $prot .= (exists $translate_p->{$codon}) ? $translate_p->{$codon} : 'X';
    }
    $prot =~ s/U$//;
    return $prot;

}

#check whether a sequence accord with gene model
#############################################
sub check_CDS{
    my $seq=shift;
    $seq =~ s/\s//g;
    $seq = uc($seq);

    my ($start,$end,$mid,$triple) = (0,0,0,0);
    my $len=length($seq);

    $triple=1 if($len%3 == 0);
    $start=1 if($seq=~/^ATG/);
    
    $end=1 if($seq=~/TAA$/ || $seq=~/TAG$/ || $seq=~/TGA$/);

    $mid=1;
    for (my $i=3; $i<$len-3; $i+=3) {
        my $codon=substr($seq,$i,3);
        $mid=0 if($codon eq 'TGA' || $codon eq 'TAG' || $codon eq 'TAA');
        print $i."\n" if($codon eq 'TGA' || $codon eq 'TAG' || $codon eq 'TAA');    
    }


    return ($start,$end,$mid,$triple);
