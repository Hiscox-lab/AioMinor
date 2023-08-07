#!/usr/bin/perl -w
use List::Util qw(sum);
$c2p{$_} = "Z" for qw(SAA);
$c2p{$_} = "J" for qw(MTT MTA);
$c2p{$_} = "B" for qw(AAT AAC GAC GAT RAT RAC);
$c2p{$_} = "L" for qw(CTA CTT CTG CTC CTN CTK TTA TTG TTR YTA YTG);
$c2p{$_} = "R" for qw(CGA CGT CGC CGG CGN AGG AGA AGR MGA MGG);
$c2p{$_} = "S" for qw(TCA TCG TCT TCC TCN AGT AGC AGY);
$c2p{$_} = "A" for qw(GCC GCT GCA GCG GCN);
$c2p{$_} = "G" for qw(GGC GGT GGA GGG GGN);
$c2p{$_} = "P" for qw(CCA CCT CCG CCC CCN);
$c2p{$_} = "T" for qw(ACA ACG ACC ACT ACN);
$c2p{$_} = "V" for qw(GTA GTC GTG GTT GTN);
$c2p{$_} = "I" for qw(ATT ATC ATY ATA ATW);
#$c2p{$_} = "_" for qw(TAA TAG TAR TGA);
$c2p{$_} = "*" for qw(TAA TAG TAR TGA);
$c2p{$_} = "C" for qw(TGT TGC TGY);
$c2p{$_} = "D" for qw(GAT GAC GAY);
$c2p{$_} = "E" for qw(GAA GAG GAR);
$c2p{$_} = "F" for qw(TTT TTC TTY);
$c2p{$_} = "H" for qw(CAT CAC CAY);
$c2p{$_} = "K" for qw(AAA AAG AAR);
$c2p{$_} = "N" for qw(AAT AAC AAY);
$c2p{$_} = "Q" for qw(CAA CAG CAR);
$c2p{$_} = "Y" for qw(TAT TAC TAY);
$c2p{$_} = "M" for qw(ATG);
$c2p{$_} = "W" for qw(TGG);
$c2p{$_} = "X" for qw(... AAN AGN ANA ANC ANG ANN ANT ATN CAN CNA CNC CNG CNN CNT GAN GNA GNC GNG GNN GNT NAA NAC NAG NAN NAT NCA NCC NCG NCN NCT NGA NGC NGG NGN NGT NNA NNC NNG NNT NTA NTC NTG NTN NTT TAN TGN TNA TNC TNG TNN TNT TTN NNN);#TYA CAM

$inputENTROPY="$ARGV[0]";
$inputREGION="$ARGV[1]";
$inputENTROPYFLITER="$ARGV[2]";
$inputCONDONALL="$ARGV[3]";
$inputAMINOACIDALL="$ARGV[4]";

$outputCONDONALLR="$ARGV[5]";
$outputAMINOACIDALLR="$ARGV[6]";
$outputMINORCHANGE="$ARGV[7]";

open(ENTROPY, "$inputENTROPY");
@entropylist=<ENTROPY>;
close ENTROPY;

foreach (@entropylist) {
    chomp;
    @eachentropy=split(/\t/);
    $hashentropyAcnt{"$eachentropy[1]\t$eachentropy[2]\t$eachentropy[3]"}=$eachentropy[6];
    $hashentropyCcnt{"$eachentropy[1]\t$eachentropy[2]\t$eachentropy[3]"}=$eachentropy[8];
    $hashentropyTcnt{"$eachentropy[1]\t$eachentropy[2]\t$eachentropy[3]"}=$eachentropy[10];
    $hashentropyGcnt{"$eachentropy[1]\t$eachentropy[2]\t$eachentropy[3]"}=$eachentropy[12];
}

open(REGION, "$inputREGION");
@regionlist=<REGION>;
close REGION;
foreach (@regionlist) {
    chomp;
    unless (/^Protein\tBeg/){
       @regions=split(/\t/);
       @begend=("$regions[1]","$regions[2]","$regions[0]");
       push (@{$hashregions{"$regions[3]"}},[@begend]);
    }
}

open(ENTROPYFLITER, "$inputENTROPYFLITER");
@entropyfilterlist=<ENTROPYFLITER>;
close ENTROPYFLITER;

foreach (@entropyfilterlist) {
    chomp;
    my @eachentropyfilter=split(/\t/);
    if ($hashentropyAcnt{"$eachentropyfilter[1]\t$eachentropyfilter[2]\t$eachentropyfilter[3]"} ne $eachentropyfilter[6]) {
        my @getApostion=split(/\t/);
        foreach my $eachprotein(@{$hashregions{$getApostion[1]}}) {
            if (@$eachprotein[0]<=$getApostion[2] && @$eachprotein[1] >= $getApostion[2]) {
                print "A\:\t$getApostion[1]\t@$eachprotein[2]\t@$eachprotein[0]\t@$eachprotein[1]\t$getApostion[2]\t";
                my $getAposincondon= ($getApostion[2] - @$eachprotein[0] + 1)/3;
                if ($getAposincondon =~/\.333/){
                   $getAposincondonvalue=int($getAposincondon)+1;
                   print $getAposincondonvalue, "\tthis is 1\n";
                   push (@collectionA,"$getApostion[1]\t@$eachprotein[2]\t$getAposincondonvalue\t0\n");
                }elsif($getAposincondon =~/\.666/){
                   $getAposincondonvalue=int($getAposincondon)+1;
                   print int($getAposincondon)+1, "\tthis is 2\n";
                   push (@collectionA,"$getApostion[1]\t@$eachprotein[2]\t$getAposincondonvalue\t1\n");
                }else{
                   print $getAposincondon, "\tthis is 3\n";
                   push (@collectionA,"$getApostion[1]\t@$eachprotein[2]\t$getAposincondon\t2\n");
                }
            } 
        }
    }
    if ($hashentropyCcnt{"$eachentropyfilter[1]\t$eachentropyfilter[2]\t$eachentropyfilter[3]"} ne $eachentropyfilter[8]) {
        my @getCpostion=split(/\t/);
        foreach my $eachprotein(@{$hashregions{$getCpostion[1]}}) {
            if (@$eachprotein[0]<=$getCpostion[2] && @$eachprotein[1] >= $getCpostion[2]) {
                print "C\:\t$getCpostion[1]\t@$eachprotein[2]\t@$eachprotein[0]\t@$eachprotein[1]\t$getCpostion[2]\t";
                my $getCposincondon= ($getCpostion[2] - @$eachprotein[0] + 1)/3;
                if ($getCposincondon =~/\.333/){
                   $getCposincondonvalue=int($getCposincondon)+1;
                   print $getCposincondonvalue, "\tthis is 1\n";
                   push (@collectionC,"$getCpostion[1]\t@$eachprotein[2]\t$getCposincondonvalue\t0\n");
                }elsif($getCposincondon =~/\.666/){
                   $getCposincondonvalue=int($getCposincondon)+1;
                   print int($getCposincondon)+1, "\tthis is 2\n";
                   push (@collectionC,"$getCpostion[1]\t@$eachprotein[2]\t$getCposincondonvalue\t1\n");
                }else{
                   print $getCposincondon, "\tthis is 3\n";
                   push (@collectionC,"$getCpostion[1]\t@$eachprotein[2]\t$getCposincondon\t2\n");
                }
            } 
        }
    }
    if ($hashentropyTcnt{"$eachentropyfilter[1]\t$eachentropyfilter[2]\t$eachentropyfilter[3]"} ne $eachentropyfilter[10]) {
        my @getTpostion=split(/\t/);
        foreach my $eachprotein(@{$hashregions{$getTpostion[1]}}) {
            if (@$eachprotein[0]<=$getTpostion[2] && @$eachprotein[1] >= $getTpostion[2]) {
                print "T\:\t$getTpostion[1]\t@$eachprotein[2]\t@$eachprotein[0]\t@$eachprotein[1]\t$getTpostion[2]\t";
                my $getTposincondon= ($getTpostion[2] - @$eachprotein[0] + 1)/3;
                if ($getTposincondon =~/\.333/){
                   $getTposincondonvalue=int($getTposincondon)+1;
                   print $getTposincondonvalue, "\tthis is 1\n";
                   push (@collectionT,"$getTpostion[1]\t@$eachprotein[2]\t$getTposincondonvalue\t0\n");
                }elsif($getTposincondon =~/\.666/){
                   $getTposincondonvalue=int($getTposincondon)+1;
                   print int($getTposincondon)+1, "\tthis is 2\n";
                   push (@collectionT,"$getTpostion[1]\t@$eachprotein[2]\t$getTposincondonvalue\t1\n");
                }else{
                   print $getTposincondon, "\tthis is 3\n";
                   push (@collectionT,"$getTpostion[1]\t@$eachprotein[2]\t$getTposincondon\t2\n");
                }
            } 
        }
    }
    if ($hashentropyGcnt{"$eachentropyfilter[1]\t$eachentropyfilter[2]\t$eachentropyfilter[3]"} ne $eachentropyfilter[12]) {
        my @getGpostion=split(/\t/);
        foreach my $eachprotein(@{$hashregions{$getGpostion[1]}}) {
            if (@$eachprotein[0]<=$getGpostion[2] && @$eachprotein[1] >= $getGpostion[2]) {
                print "G\:\t$getGpostion[1]\t@$eachprotein[2]\t@$eachprotein[0]\t@$eachprotein[1]\t$getGpostion[2]\t";
                my $getGposincondon= ($getGpostion[2] - @$eachprotein[0] + 1)/3;
                if ($getGposincondon =~/\.333/){
                   $getGposincondonvalue=int($getGposincondon)+1;
                   print $getGposincondonvalue, "\tthis is 1\n";
                   push (@collectionG,"$getGpostion[1]\t@$eachprotein[2]\t$getGposincondonvalue\t0\n");
                }elsif($getGposincondon =~/\.666/){
                   $getGposincondonvalue=int($getGposincondon)+1;
                   print int($getGposincondon)+1, "\tthis is 2\n";
                   push (@collectionG,"$getGpostion[1]\t@$eachprotein[2]\t$getGposincondonvalue\t1\n");
                }else{
                   print $getGposincondon, "\tthis is 3\n";
                   push (@collectionG,"$getGpostion[1]\t@$eachprotein[2]\t$getGposincondon\t2\n");
                }
            } 
        }
    }    
}


open(CONDONALL, "$inputCONDONALL");
@condonalllist=<CONDONALL>;
@condonalllistr=@condonalllist;
close CONDONALL;

foreach (@condonalllistr) {
    chomp;
    my @eachcondonall=split(/\t/);
    my @topcodon=split(//,$eachcondonall[7]);
    push (@topcodon,"$eachcondonall[4]");
    push (@topcodon,"$eachcondonall[5]");
    push (@topcodon,"$eachcondonall[6]");
    push (@topcodon,"$eachcondonall[7]");
    push (@topcodon,"$eachcondonall[8]");
    push (@topcodon,"$eachcondonall[9]");
    push (@{$hashcondonall{"$eachcondonall[0]\t$eachcondonall[1]\t$eachcondonall[2]"}},[@topcodon]);
}

open(AMINOACIDALL, "$inputAMINOACIDALL");
@aminoacidlist=<AMINOACIDALL>;
close AMINOACIDALL;

foreach (@collectionA) {
    chomp;
    my @splitcollectionA=split(/\t/);
    @aminoacidlist= grep (!/^$splitcollectionA[0]\t$splitcollectionA[1]\t$splitcollectionA[2]\t/, @aminoacidlist);
    @condonalllist= grep (!/^$splitcollectionA[0]\t$splitcollectionA[1]\t$splitcollectionA[2]\t/, @condonalllist);
    my @foundlist = grep(grep(!/A/, @$_[$splitcollectionA[3]]), @{$hashcondonall{"$splitcollectionA[0]\t$splitcollectionA[1]\t$splitcollectionA[2]"}});
    my $num=0;
    my %hashtogetnewaanum; my %hashtogetnewaanumnew;
    foreach (@foundlist) {
        $num++;
        #print "$splitcollectionA[0]\t$splitcollectionA[1]\t$splitcollectionA[2]\t", @$_, "\n";
        my $raa= $c2p{uc("@$_[0]@$_[1]@$_[2]")};
        push (@condonalllist, "$splitcollectionA[0]\t$splitcollectionA[1]\t$splitcollectionA[2]\t$num\t@$_[3]\t@$_[4]\t@$_[5]\t@$_[6]\t@$_[7]\t@$_[8]\n");
        push (@{$hashtogetnewaanum{"$raa"}},@$_[7]);
        $refaain=@$_[3];
        $refsitein=@$_[4];
        $refcodonin=@$_[5];
        $aacoveragein=@$_[8];
    }
    
    for my $key (keys %hashtogetnewaanum) {
        #print "$key => ";
        my $value_arr = $hashtogetnewaanum{$key};
        $hashtogetnewaanumnew{$key}=sum(@$value_arr);
    }
    
    my $nump=0;
    foreach $newkyes (sort {$hashtogetnewaanumnew{$b} <=> $hashtogetnewaanumnew{$a}} keys %hashtogetnewaanumnew) {
        $nump++;
        push (@aminoacidlist, "$splitcollectionA[0]\t$splitcollectionA[1]\t$splitcollectionA[2]\t$nump\t$refaain\t$refsitein\t$refcodonin\t$newkyes\t$hashtogetnewaanumnew{$newkyes}\t$aacoveragein\n");
    }
}

foreach (@collectionC) {
    chomp;
    my @splitcollectionC=split(/\t/);
    @aminoacidlist= grep (!/^$splitcollectionC[0]\t$splitcollectionC[1]\t$splitcollectionC[2]\t/, @aminoacidlist);
    @condonalllist= grep (!/^$splitcollectionC[0]\t$splitcollectionC[1]\t$splitcollectionC[2]\t/, @condonalllist);
    my @foundlist = grep(grep(!/C/, @$_[$splitcollectionC[3]]), @{$hashcondonall{"$splitcollectionC[0]\t$splitcollectionC[1]\t$splitcollectionC[2]"}});
    my $num=0;
    my %hashtogetnewaanum; my %hashtogetnewaanumnew;
    foreach (@foundlist) {
        $num++;
        #print "$splitcollectionC[0]\t$splitcollectionC[1]\t$splitcollectionC[2]\t", @$_, "\n";
        my $raa= $c2p{uc("@$_[0]@$_[1]@$_[2]")};
        push (@condonalllist, "$splitcollectionC[0]\t$splitcollectionC[1]\t$splitcollectionC[2]\t$num\t@$_[3]\t@$_[4]\t@$_[5]\t@$_[6]\t@$_[7]\t@$_[8]\n");
        push (@{$hashtogetnewaanum{"$raa"}},@$_[7]);
        $refaain=@$_[3];
        $refsitein=@$_[4];
        $refcodonin=@$_[5];
        $aacoveragein=@$_[8];
    }
    
    for my $key (keys %hashtogetnewaanum) {
        #print "$key => ";
        my $value_arr = $hashtogetnewaanum{$key};
        $hashtogetnewaanumnew{$key}=sum(@$value_arr);
    }
    
    my $nump=0;
    foreach $newkyes (sort {$hashtogetnewaanumnew{$b} <=> $hashtogetnewaanumnew{$a}} keys %hashtogetnewaanumnew) {
        $nump++;
        push (@aminoacidlist, "$splitcollectionC[0]\t$splitcollectionC[1]\t$splitcollectionC[2]\t$nump\t$refaain\t$refsitein\t$refcodonin\t$newkyes\t$hashtogetnewaanumnew{$newkyes}\t$aacoveragein\n");
    }
}

foreach (@collectionT) {
    chomp;
    my @splitcollectionT=split(/\t/);
    @aminoacidlist= grep (!/^$splitcollectionT[0]\t$splitcollectionT[1]\t$splitcollectionT[2]\t/, @aminoacidlist);
    @condonalllist= grep (!/^$splitcollectionT[0]\t$splitcollectionT[1]\t$splitcollectionT[2]\t/, @condonalllist);
    my @foundlist = grep(grep(!/T/, @$_[$splitcollectionT[3]]), @{$hashcondonall{"$splitcollectionT[0]\t$splitcollectionT[1]\t$splitcollectionT[2]"}});
    my $num=0;
    my %hashtogetnewaanum; my %hashtogetnewaanumnew;
    foreach (@foundlist) {
        $num++;
        #print "$splitcollectionT[0]\t$splitcollectionT[1]\t$splitcollectionT[2]\t", @$_, "\n";
        my $raa= $c2p{uc("@$_[0]@$_[1]@$_[2]")};
        push (@condonalllist, "$splitcollectionT[0]\t$splitcollectionT[1]\t$splitcollectionT[2]\t$num\t@$_[3]\t@$_[4]\t@$_[5]\t@$_[6]\t@$_[7]\t@$_[8]\n");
        push (@{$hashtogetnewaanum{"$raa"}},@$_[7]);
        $refaain=@$_[3];
        $refsitein=@$_[4];
        $refcodonin=@$_[5];
        $aacoveragein=@$_[8];
    }
    
    for my $key (keys %hashtogetnewaanum) {
        #print "$key => ";
        my $value_arr = $hashtogetnewaanum{$key};
        $hashtogetnewaanumnew{$key}=sum(@$value_arr);
    }
    
    my $nump=0;
    foreach $newkyes (sort {$hashtogetnewaanumnew{$b} <=> $hashtogetnewaanumnew{$a}} keys %hashtogetnewaanumnew) {
        $nump++;
        push (@aminoacidlist, "$splitcollectionT[0]\t$splitcollectionT[1]\t$splitcollectionT[2]\t$nump\t$refaain\t$refsitein\t$refcodonin\t$newkyes\t$hashtogetnewaanumnew{$newkyes}\t$aacoveragein\n");
    }
}

foreach (@collectionG) {
    chomp;
    my @splitcollectionG=split(/\t/);
    @aminoacidlist= grep (!/^$splitcollectionG[0]\t$splitcollectionG[1]\t$splitcollectionG[2]\t/, @aminoacidlist);
    @condonalllist= grep (!/^$splitcollectionG[0]\t$splitcollectionG[1]\t$splitcollectionG[2]\t/, @condonalllist);
    my @foundlist = grep(grep(!/G/, @$_[$splitcollectionG[3]]), @{$hashcondonall{"$splitcollectionG[0]\t$splitcollectionG[1]\t$splitcollectionG[2]"}});
    my $num=0;
    my %hashtogetnewaanum; my %hashtogetnewaanumnew;
    foreach (@foundlist) {
        $num++;
        #print "$splitcollectionG[0]\t$splitcollectionG[1]\t$splitcollectionG[2]\t", @$_, "\n";
        my $raa= $c2p{uc("@$_[0]@$_[1]@$_[2]")};
        push (@condonalllist, "$splitcollectionG[0]\t$splitcollectionG[1]\t$splitcollectionG[2]\t$num\t@$_[3]\t@$_[4]\t@$_[5]\t@$_[6]\t@$_[7]\t@$_[8]\n");
        push (@{$hashtogetnewaanum{"$raa"}},@$_[7]);
        $refaain=@$_[3];
        $refsitein=@$_[4];
        $refcodonin=@$_[5];
        $aacoveragein=@$_[8];
    }
    
    for my $key (keys %hashtogetnewaanum) {
        #print "$key => ";
        my $value_arr = $hashtogetnewaanum{$key};
        $hashtogetnewaanumnew{$key}=sum(@$value_arr);
    }
    
    my $nump=0;
    foreach $newkyes (sort {$hashtogetnewaanumnew{$b} <=> $hashtogetnewaanumnew{$a}} keys %hashtogetnewaanumnew) {
        $nump++;
        push (@aminoacidlist, "$splitcollectionG[0]\t$splitcollectionG[1]\t$splitcollectionG[2]\t$nump\t$refaain\t$refsitein\t$refcodonin\t$newkyes\t$hashtogetnewaanumnew{$newkyes}\t$aacoveragein\n");
    }
}

open(CONDONALLR, ">$outputCONDONALLR");
open(AMINOACIDALLR, ">$outputAMINOACIDALLR");
foreach (@condonalllist) {
    chomp;
    if (/^Chr\tProtein/) {
        print CONDONALLR $_,"\n";
    }else{
        my @condontables=split(/\t/);
        my @groups;
        push (@groups, "$condontables[2]");
        push (@groups, "$condontables[3]");
        $hashcondontables{"$condontables[0]\t$condontables[1]\t$condontables[2]\t$condontables[3]\t$condontables[4]\t$condontables[5]\t$condontables[6]\t$condontables[7]\t$condontables[8]\t$condontables[9]"}=[@groups];
        my @groupelecments=("$condontables[1]","$condontables[2]","$condontables[3]","$condontables[4]","$condontables[5]","$condontables[6]","$condontables[7]","$condontables[8]","$condontables[9]");
        push (@{$hashelecments{"$condontables[1]\t$condontables[2]\t$condontables[4]\t$condontables[9]"}}, [@groupelecments]);
    }
}

foreach (@aminoacidlist) {
    chomp;
    if (/^Chr\tProtein/) {
        print AMINOACIDALLR $_,"\n";
    }else{
        my @aminoacidtables=split(/\t/);
        my @groups;
        push (@groups, "$aminoacidtables[2]");
        push (@groups, "$aminoacidtables[3]");
        $hashaminoacidtable{"$aminoacidtables[0]\t$aminoacidtables[1]\t$aminoacidtables[2]\t$aminoacidtables[3]\t$aminoacidtables[4]\t$aminoacidtables[5]\t$aminoacidtables[6]\t$aminoacidtables[7]\t$aminoacidtables[8]\t$aminoacidtables[9]"}=[@groups];
    }
}

open(MINORCHANGE, ">$outputMINORCHANGE");
print MINORCHANGE "Protein\tAAPosition\tRefAA\tAAcoverage\tCntNonSyn\tCntSyn\tCntStop\tCntN\n";
shift @regionlist;
foreach (@regionlist) {
    chomp;
    @eachcodongregions=split(/\t/);
    foreach my $newkyes (sort {$hashcondontables{$a}[0] <=> $hashcondontables{$b}[0] or $hashcondontables{$a}[1] <=> $hashcondontables{$b}[1]} keys %hashcondontables) {
        if ($newkyes=~/^$eachcodongregions[3]\t$eachcodongregions[0]\t/) {
            print CONDONALLR $newkyes,"\n";
        }
    }
                
    foreach my $newkyes (sort {$hashaminoacidtable{$a}[0] <=> $hashaminoacidtable{$b}[0] or $hashaminoacidtable{$a}[1] <=> $hashaminoacidtable{$b}[1]} keys %hashaminoacidtable) {
        if ($newkyes=~/^$eachcodongregions[3]\t$eachcodongregions[0]\t/) {
            print AMINOACIDALLR $newkyes,"\n";
        }
    }

    foreach my $key1 (sort {$hashelecments{$a}[0][1] <=> $hashelecments{$b}[0][1]} keys %hashelecments) {
        if ($key1=~/^$eachcodongregions[0]\t/) {
            #print "$key1 => ";
            print MINORCHANGE "$key1\t";
            @aaref=split(/\t/, $key1);
            my @codoncollection;
            my %hashaacount;
            foreach my $elementingroup(@{$hashelecments{"$key1"}}) {
                unless (@$elementingroup[5] eq @$elementingroup[6]) {
                    #print "@$elementingroup[5]\t@$elementingroup[6]";
                    push (@{$hashaacount{$c2p{uc(@$elementingroup[6])}}}, "@$elementingroup[7]");
                    #push (@codoncollection, $hashaacount{$c2p{uc(@$elementingroup[6])}});
                }
            }
            #print "@codoncollection\n";
            my @nonsym;
            my @sym;
            my @countN;
            my @stopnum;
            for my $key2 (keys %hashaacount) {
                #print "$key2->";
                my $value_arr = $hashaacount{$key2};
                #print join(",", @$value_arr), "\ ";
                if ($key2 eq $aaref[2]) {
                    push (@sym, sum(@$value_arr));
                }
                if ($key2 eq "X") {
                    push (@countN, sum(@$value_arr));
                }
                if ($key2 ne $aaref[2] && $key2 ne "X") {
                    push (@nonsym, sum(@$value_arr));
                }
                if ($key2 eq "\*") {
                    push (@stopnum, sum(@$value_arr));
                }
            }
    
            if ($#nonsym >= 0) {
                #print "nonsym\:", sum(@nonsym),"\ ";
                print MINORCHANGE sum(@nonsym),"\t";
            }else{
                #print "nonsym\:", 0 ,"\ ";
                print MINORCHANGE "0\t";
            }

            if ($#sym >= 0) {
                #print "\:\:\: sym\:", sum(@sym),"\ ";
                print MINORCHANGE sum(@sym),"\t";
            }else{
                #print "\:\:\: sym\:", 0 ,"\ ";
                print MINORCHANGE "0\t";
            }

            if ($#stopnum >= 0) {
                #print "Stop\:", sum(@stopnum),"\ ";
                print MINORCHANGE sum(@stopnum),"\t";
            }else{
                #print "Stop\:", 0 ,"\ ";
                print MINORCHANGE "0\t";
            }

            if ($#countN >= 0) {
                #print "N\:", sum(@countN),"\ ";
                print MINORCHANGE sum(@countN),"\t";
            }else{
                #print "N\:", 0 ,"\ ";
                print MINORCHANGE "0\t";
            }
            #print "\n";
            print MINORCHANGE "\n";
        }
    }
}
close CONDONALLR;
close AMINOACIDALLR;
