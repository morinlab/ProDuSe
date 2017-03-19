#!/usr/bin/perl
use strict;
#stream the raw Produse VCF through this script to filter on dual strand support and strand bias features
my $p_thresh = 0.05;
my %ixmap = ("A"=>0,"C"=>1,"G"=>2,"T"=>3);

my $vaf_thresh = $ARGV[0] || 0.01; #set to default VAF cutoff of 1% if no threshold specified

my $min_duplex = 3;

#if not duplex
my $min_SP = 1;
my $min_SN = 1;
my $min_tot = 3;

while(my $line = <STDIN>){
    chomp $line;
    my $skip;
    if ($line =~ /#/) {
        print "$line\n";
        next;
    }
    my @f = split /\t/, $line;
    my $ref = $f[3];
    my $nref = $f[4]; #may be a list of SNVs. Skip all of these?
    next if $nref =~ /,/;
    my $ref_ix = $ixmap{$ref};
    my $nref_ix = $ixmap{$nref};
    my @details = split /;/, $f[7];
    my $duplex_PN_count = 0;
    my $SP_count;
    my $SN_count;
    my $vaf;
    for my $set (@details){
	if ($set =~ /StrBiasP=(.+)/){
	    my @ps = split /,/, $1;
	    my $p = $ps[$nref_ix];
	    if($p < $p_thresh){
	    $skip++;
	    }
	}
	if($set =~ /DPN=(.+)/){
	    my @ps = split /,/, $1;
	    my $p = $ps[$nref_ix];
	    $duplex_PN_count = $p;
	}
	
	if($set =~ /SP=(.+)/){
	    my @ps = split /,/, $1;
	    my $p = $ps[$nref_ix];
	    $SP_count = $p;
	}
	if($set =~ /SN=(.+)/){
	    my @ps = split /,/, $1;
	    my $p = $ps[$nref_ix];
	    $SN_count =$p;
	}
	if($set =~ /MC=(.+)/){
	    my @ps = split /,/, $1;
	    my $mut = $ps[$nref_ix];
	    my $wt = $ps[$ref_ix];
	    $vaf = $mut/($wt + $mut);
	    if($vaf < $vaf_thresh){
		$skip = 1;
	    }
	}
	
    }
    

    next if $skip;
    if($duplex_PN_count < $min_duplex){
	unless($SP_count >= $min_SP && $SN_count >= $min_SN){
	    next;
	}
	unless ($SP_count + $SN_count >= $min_tot){
	    next;
	}
    }
    
    print "$line,VAF=$vaf\n";
}
