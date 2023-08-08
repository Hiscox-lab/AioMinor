#!/usr/bin/perl
# Calculates the size of fasta files
@files = @ARGV;
if ($files[1] eq "fix"){
	print "\nOK, fixing sam file.\n";
	&fix_mac($files[0]);
}
$POS = 0;
$genome = "";
$test = -1;
$G = 0;
$A = 0;
$T = 0;
$C = 0;
open(INFILE, "$ARGV[0]"); # opens quasirecomd output
open(INCOVERAGE, "$ARGV[1]"); # opens quasirecomd output
@eachcoverage=<INCOVERAGE>;
@numposcov=split(/\t/,$eachcoverage[0]);
close INCOVERAGE;

open(OUTB, ">$ARGV[2]");
print OUTB ">consensus\n";
print OUTB "N" x $numposcov[0];

while($line = <INFILE>){
chomp $line;
unless ($line=~/#Offset|^Pos/) {
	# split the line into an array called cells
	($POS,$A,$C,$G,$T,$gap) = split(/\t/, $line);
    $A = $A * 1; $G = $G * 1; $C = $C * 1; $T = $T * 1;
    #print "$POS\t$a\t$A\t$C\t$G\t$T\n";
    $test = -1;
    $consensus = "N";
    if ($A > $test){$test = $A; $consensus = "A";}
    if ($C > $test){$test = $C; $consensus = "C";}
    if ($G > $test){$test = $G; $consensus = "G";}
    if ($T > $test){$test = $T; $consensus = "T";}
    $genome = $genome.$consensus;
    #print "$POS\t$consensus\n";
    }
}
print OUTB "$genome\n";
	
close(OUTB);

sub fix_mac
{
	$filefix = $_[0];
	chomp $filefix;
	print "Correcting a file.\n";
	open (OUT, ">temp");
	open (INFILE, "$filefix");
	while (<INFILE>){s/[\r\n]+/\n/g; print OUT "$_";}
	close OUT;
	close INFILE;
	system("mv temp $filefix"); 
}
