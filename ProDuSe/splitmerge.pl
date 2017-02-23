#!/usr/bin/perl
use strict;
my $n = 0;
my $m = 0;
my $cap = chr(69+33);
while(<STDIN>){
    chomp;
    my @fields = split /\t/;

    my $tag = $fields[11];
    #e.g. 17F118S17R
    unless($tag =~ /XD/){
	#no tricks needed on positions. Just update flags to be consistent (can this be done?)
	print "$_\n";
	next;
    }
    my ($f,$b,$r);
    unless($tag =~ /R/){
	#use forward
	$n = 1;
	$r  = 0;
	$b = length($fields[9]);
	$f = 0;
	$tag = "";
#	print "$n $r $b $f\n";
    }
    unless($tag =~ /F/){
	#use reverse
	$n = 0;
	$r  = 0;
	$b = length($fields[9]);
	$f = 0;
	$tag = "";
#	print "$n $r $b $f\n";
    }
    if($tag =~ /:\d+S/ && $tag =~ /\d+S$/ ){
#	print "MATCH HERE\n";
	    $b = length($fields[9]);
	    $f = 0;
	    $r = 0;
	    $m++;
	    $tag = "";
    }
    if($tag =~/\d+F\d+S\d+F\d+S/){
#	print "M2\n";
	while($tag =~ /(\d+)[F]/g){
	    $f+=$1;
	    #print STDERR "F:$1, $f\n";
	}
	while($tag =~ /(\d+)[S]/g){
	    $b+=$1;
	    #print STDERR "B:$1, $b\n";
	}
	while($tag =~ /(\d+)[R]/g){
	    $r+=$1;
	}
	my $len = length($fields[9]);
#	print STDERR "$len - $b - $f $fields[9]<\n";
	$r = $len - $b - $f;
	$n = 1; #set back to odd value
	#print STDERR "TAG: $tag\n$f\n$b\n$r (give to F)\n";
    }
    elsif($tag =~/\d+F\d+S\d+R\d+S/){
#	print "M3\n";
	while($tag =~ /(\d+)[F]/g){
	    $f+=$1;
	}
	while($tag =~ /(\d+)[S]/g){
	    $b+=$1;
	}
	while($tag =~ /(\d+)[R]/g){
	    $r+=$1;
	}
	my $len= length($fields[9]);
	$b = $len - $r - $f;
	$n = 2; #set back to even value
	#print STDERR "TAG: $tag\n$f\n$b\n$r (give to R)\n";
    }
    elsif($tag =~ /(\d+)F(\d+)S(\d+)R/){
	#next;
#	print "M4\n";
	($f,$b,$r) = ($1,$2,$3);
	#print "$f $b $r\n";
    }
    elsif($tag =~ /(\d+)S(\d+)R/){
	($f,$b,$r) = (0,$1,$2);
#	print "M5\n";
	#print "$f $b $r\n";
	#next;
    }
    elsif($tag =~ /(\d+)F(\d+)S/){
	($f,$b,$r) = ($1,$2,0);
#	print "M6\n";
	#print "$f $b $r\n";
	#next;
    }
    elsif($tag =~ /(\d+)S/){
#	print "M7\n";
	($f,$b,$r) = (0,$1,0);
	#print "$f $b $r\n";
	#next;
    }
 
    #randomly (but consistently) take the first + shared or last + shared region as read 1 or 2 and the remainder for the other
    
    my $seq = $fields[9];
    my $qual = $fields[10];
    my $cigar = $fields[5];
    #next unless $cigar =~/I/;
    my ($l1,$l2);
    if($n % 2){ #odd, give to forward
	$l1 = $f + $b;
	$l2 = $r;
	#print "$l1 and $l2\n";
    }
    else{
	
	$l1 = $f;
	$l2 = $r + $b;
    }
    
 #   print "$tag: $l1 and $l2 ($f,$b,$r)\n";
    my $seq1 = substr $seq, 0, $l1;
    my $seq2 = substr $seq, -$l2 ;
  #  print "S1:$seq1\nS2:$seq2\n$f and $r\n";
    #if(length($seq2) < $l2){
#	print "$seq2 < $l2\n";
#    }
    my $qual1 = substr $qual, 0, $l1;
    my $qual2 = substr $qual, -$l2;
    my $pos1 = $fields[3];
    
    my $flag1;
    my $flag2;
    if($fields[0]=~ /:-:/){
	#minus strand PCR family
	$flag1 = 163;
	$flag2 = 83;
	
    }
    else{
	$flag1 = 99;
	$flag2 = 147;
	
    }
    my $cig1;
    my $cig1l;
    my $cig2;
    my $cig2l;
    my $donec1;
    my $skip;
    my $genome_base_offset = 0; #shifts by any M and I operation (but not D) in the read 1
    if($f == 0 && $r == 0){
	$genome_base_offset = 0;
	#just recycle the same cigar and move on
	if($m % 2){
	    #forward read only
	    $cig2l = 0;
	    $cig1 = $cigar;
	    $cig1l = 100;
	    $seq1 = $seq;
	    $qual1 = $qual;
	}
	else{
	    $cig1l = 0;
	    $cig2l = 100;
	    $cig2 = $cigar;
	    $seq2 = $seq;
	    $qual2 = $qual;
	}
    }else{
	
	while($cigar =~ /(\d+)([MISD])/g){
	    
	my $op = $2;
	my $bases = $1;
	if($op eq "S"){
	    $skip = 1;
	}
	if($op eq "I"){
	    $skip = 1;
	}
	if($op eq "D"){
	    $skip = 1;
	    #no loop needed because these are not bases in the read
	    if($donec1){
		$cig2 .= "$bases$op"
	    }
	    else{
		$cig1 .= "$bases$op";
	    }
	    next;
	}
	#print "OP: $op, B: $bases\n";
	my $thisL;
	if($op eq "I"){
	    $genome_base_offset +=$bases;	    
	}
	while($bases>0){
		unless($donec1){
		    $genome_base_offset++;
		}
	    
	    if($cig1l<$l1){
		$cig1l++;
		$bases--;
		
		$thisL++;
		#print "cig1l $cig1l cig1 $cig1 thisL $thisL op $op\n";
		
	    }
	    elsif($cig1l==$l1 and !$donec1){
		#print "cig1l $cig1l cig1 $cig1 thisL $thisL op $op\n";
		$donec1 = 1;
		    
		$cig1 .= $thisL;
		$cig1 .= $op;

		#add these to the cigar string to complete it and start the next one
		$thisL = 1;
		
		
		$cig2l++;
		$bases--;
		#print "cig1l $cig1l cig1 $cig1 thisL $thisL op $op $cig1\n";
		#print "cig2l $cig2l cig2 $cig2 thisL $thisL op $op\n";
	    }
	    else{
		$cig2l++;
		$bases--;
		#print "cig2l $cig2l cig2 $cig2 thisL $thisL op $op\n";
		$thisL++;
		
		
	    }
	}
	if($donec1){
	    $cig2.=$thisL;
	    $cig2.=$op;
	    #print "cig2l $cig2l cig2 $cig2 thisL $thisL op $op\n";
	}
	else{
	    
	    $cig1.=$thisL;
	    $cig1.=$op;
	}
	}
    }
    $n++;
    if($f == 0 && $r == 0){
	
    }
    else{
	$genome_base_offset = length($seq1); #need to add number of D bases once that gets fixed
    }
    my $pos2 = $pos1 + $genome_base_offset;
    #print "C1L $cig1l\nC1 $cig1\nC2L $cig2l\nC2 $cig2\n";
    next if $skip;
    #print "$cig1l and $cig2l (cigars)\n";
    my @sam1 = @fields;
    my @sam2 = @fields;
    $sam1[1] = $flag1;
   
    $sam2[1] = $flag2;
    $sam1[3] = $pos1;
    $sam2[3] = $pos2;
    $sam1[6] = $sam2[2];
    $sam2[6] = $sam1[2];
    $sam1[7] = $pos2;
    $sam2[7] = $pos1;
    $sam1[9] = $seq1;
    $sam2[9] = $seq2;
    if($cap){
	$sam1[10] = capqual($qual1);
    }
    else{
	$sam1[10] = $qual1;
    }
    if($cap){
	$sam2[10] = capqual($qual2);
    }
    else{
	$sam2[10] = $qual2;
    }
    $sam1[5] = $cig1;
    $sam2[5] = $cig2;
    
    if($cig1l > 1 && length($seq1)>2){
	
	print join "\t", @sam1, "\n";
    }
    if($cig2l > 1 && length($seq2) > 2){
	print join "\t", @sam2, "\n";
    }
}
    
sub capqual{
    
    my $qstring = shift;
    my @quals = split //, $qstring;
    my $newqual;
    for(@quals){
	my $q = ord($_)-33;
	if ($q > 69){
	    #print "$_ too high at $q: swap for $cap\n";
	    $newqual.=$cap;
	}
	else{
	    $newqual.=$_;
	}
    }
    #print "from: $qstring\nto: $newqual\n";
    return $newqual;
}
