#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw(sum);

my $selfpath=dirname(__FILE__);

my %options;
my @standard_options =("help|h!",
                       "version|v!",
                       "platform=s",
                       "method=s",
                       "maxins=s",
                       "minins=s",
                       "ref=s",
                       "codingRegion=s",
                       "fq=s",
                       "fq1=s",
                       "fq2=s",
                       "primerbed=s",
                       "o=s",
                       "primer_bed=s",
                       "t|thread=s",
                       "samplename=s"
                       );

GetOptions( \%options, @standard_options );

###################### parameters setting ######################
# if no arguments supplied print the usage and exit
if ((keys (%options))==0) {
    print "please use -h to get help\n";
    exit;
}

# If the -help option is set, print the usage and exit
if ($options{'help'}) {
    print "\nUsage example\:
  perl AioMinor.pl -platform illumina -method cDNA -ref genome.fasta -codingRegion CodingRegion.txt -fq1 1.fastq.gz -fq2 2.fastq.gz
  perl AioMinor.pl -platform illumina -method amplicon -maxins 500 -ref genome.fasta -codingRegion CodingRegion.txt -primerbed primer.bed -fq1 1.fastq.gz -fq2 2.fastq.gz
  perl AioMinor.pl -platform nanopore -method amplicon -maxins 800 -ref genome.fasta -codingRegion CodingRegion.txt -primerbed primer.bed -fq fastq.gz
  
Required options:
  -platform           sequencing platform \"nanopore\" or \"illumina\".
  -method             method of sequencing Library Preparation \"amplicon\" or \"cDNA\".
  -maxins             maximum fragment length of in your amplicon Library with paired-end sequencing.
  -ref                reference genome sequence in fasta file.
  -codingRegion       coding regions in your reference genome.
  -primerbed          primer scheme in bed file.
  -fq                 fastq(.gz) file for single-end sequencing.
  -fq1                fastq(.gz) file for paired-end sequencing mate 1s.
  -fq2                fastq(.gz) file for paired-end sequencing mate 2s.
  -primerbed          amplicon primer information in bed file.

Optional options:
  -samplename         sample name, \"alignment\" by default.
  -t/-thread          number of threads, 1 by default.
  -minins             minimum fragment length of in your amplicon Library with paired-end sequencing, \"50\" by default.
  -o                  output path, \"./\" by default.
  
  -v/-version         Print version information.
  -h/-help            Print help message.\n\n";
    exit;
}
#The maximum fragment length for valid paired-end alignments.

if ($options{'version'}) {
    print "v0.02\n";
    exit;
}

###########start to process
if (!exists $options{'platform'}) {
    print "Please provide sequencing platform\n";
    exit;
}
if (!exists $options{'method'}) {
    print "Please provide method of sequencing Library Preparation\n";
    exit;
}
if (!exists $options{'ref'}) {
    print "Please provide reference genome sequence\n";
    exit;
}

my $outputpath;
if (!exists ($options{'o'})) {
    $outputpath= "\./";
}elsif (exists ($options{'o'})) {
    $outputpath= $options{'o'};
    mkdir "$outputpath";
}

my $thread;
if (!exists ($options{'t'})) {
    $thread= 1;
}elsif (exists ($options{'t'})) {
    $thread= $options{'t'};
}

my $list;
if (!exists ($options{'samplename'})) {
    $list= "alignment";
}elsif (exists ($options{'samplename'})) {
    $list= $options{'samplename'};
}

my $minins;
if (!exists ($options{'minins'})) {
    $minins= 50;
}elsif (exists ($options{'minins'})) {
    $minins= $options{'minins'};
}

if ($options{'platform'} eq "illumina" && $options{'method'} eq "amplicon") {
    print "The sequencing platform is illumina\n";
    print "The sequencing library is amplicon\n";
    
    if (!exists $options{'codingRegion'}) {
        print "Please provide coding regions in your reference genome";
        exit;
    }

    unless (exists $options{'fq'} or exists $options{'fq1'} or exists $options{'fq2'}) {
        print "Please provide the fastq files\n";
        exit;
    }
    
    if (!exists $options{'primerbed'}) {
        print "Please provide primer scheme in bed file\n";
        exit;
    }
    &primerbed;
    &illuminaamplicon;
}

if ($options{'platform'} eq "nanopore" && $options{'method'} eq "amplicon") {
    print "The sequencing platform is nanopore\n";
    print "The sequencing library is amplicon\n";
    
    if (!exists $options{'codingRegion'}) {
        print "Please provide coding regions in your reference genome";
        exit;
    }

    if (!exists $options{'fq'} or exists $options{'fq1'} or exists $options{'fq2'}) {
        print "Please provide a single basecalled fastq file with \"-fq\"\n";
        exit;
    }
    
    if (!exists $options{'primerbed'}) {
        print "Please provide primer scheme in bed file\n";
        exit;
    }
    &nanoporeamplicon;
}

if ($options{'platform'} eq "illumina" && $options{'method'} eq "cDNA") {
    print "The sequencing platform is illumina\n";
    print "The sequencing library is cDNA\n";
    
    if (!exists $options{'codingRegion'}) {
        print "Please provide coding regions in your reference genome";
        exit;
    }

    unless (exists $options{'fq'} or exists $options{'fq1'} or exists $options{'fq2'}) {
        print "Please provide the fastq files\n";
        exit;
    }
    
    if (exists $options{'primerbed'}) {
        print "\"-primerbed\" is only for the ampilcon sequencing\n";
        exit;
    }
    
    if (exists $options{'maxins'}) {
        print "\"-maxins\" is only for the paired-end illumina ampilcon sequencing\n";
        exit;
    }

    &illuminacdna;
}

sub illuminacdna {
    ###mapping run 1
    print "The sample name is $list\n";
    my $indexfileExist = -e "$options{'ref'}\.1.bt2";
    if ($indexfileExist) {
        print "The bowtie2-build has been done\n";
    } else {
        print "The bowtie2-build is running\n";
        system ("bowtie2-build $options{'ref'} $options{'ref'}");
    }
    
    mkdir "$outputpath/alignment";
    if (exists $options{'fq'} and !exists $options{'fq1'} and !exists $options{'fq2'}) {
        print $options{'fq'}, "\n";
        system("bowtie2 --local -x $options{'ref'} -p $thread -U $options{'fq'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2308 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }else {
        system("bowtie2 --local -X 2000 --no-mixed -x $options{'ref'} -p $thread -1 $options{'fq1'} -2 $options{'fq2'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2316 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }
    unlink ("$outputpath/alignment/$list\.sam");
    system("samtools index $outputpath/alignment/$list\.sorted.bam");
    system("java -jar $selfpath/Scripts/picard.jar MarkDuplicates INPUT=$outputpath/alignment/$list\.sorted.bam OUTPUT=$outputpath/alignment/$list\_markup.bam METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=tmp ASSUME_SORTED=true");
    system("samtools index $outputpath/alignment/$list\_markup.bam");
    unlink ("$outputpath/alignment/$list\.sorted.bam");
    unlink ("$outputpath/alignment/$list\.sorted.bam.bai");
    system("java -XX:+UseParallelGC -XX:+UseNUMA -XX:NewRatio=9 -Xms10G -Xmx200G -jar $selfpath/Scripts/QuasiRecomb.jar -i $outputpath/alignment/$list\_markup.bam -o $outputpath/alignment/$list\.support_oem -coverage");
    system("perl $selfpath/Scripts/parse_quasirecomb_V2_new.pl $outputpath/alignment/$list\.support_oem/support/allel_distribution_phred_weighted.txt $outputpath/alignment/$list\.support_oem/support/coverage.txt $outputpath/alignment/$list\.support_oem/consensus.txt");
    system ("rm -r $outputpath/alignment/$list\.support_oem/support");

    ###mapping run 2
    system("bowtie2-build --threads $thread $outputpath/alignment/$list\.support_oem/consensus.txt $outputpath/alignment/$list\.support_oem/consensus");
    if (exists $options{'fq'} and !exists $options{'fq1'} and !exists $options{'fq2'}) {
        print $options{'fq'}, "\n";
        system("bowtie2 --local -x $outputpath/alignment/$list\.support_oem/consensus -p $thread -U $options{'fq'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2308 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }else {
        system("bowtie2 --local -X 2000 --no-mixed -x $outputpath/alignment/$list\.support_oem/consensus -p $thread -1 $options{'fq1'} -2 $options{'fq2'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2316 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }
    unlink ("$outputpath/alignment/$list\.sam");
    system("samtools index $outputpath/alignment/$list\.sorted.bam");
    system("java -jar $selfpath/Scripts/picard.jar MarkDuplicates INPUT=$outputpath/alignment/$list\.sorted.bam OUTPUT=$outputpath/alignment/$list\_markup.bam METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=tmp ASSUME_SORTED=true");
    system("samtools index $outputpath/alignment/$list\_markup.bam");
    unlink ("$outputpath/alignment/$list\.sorted.bam");
    unlink ("$outputpath/alignment/$list\.sorted.bam.bai");

    ###running diversiutils, filtration, minor variation
    mkdir "$outputpath/1_Syn_NonSyn_aa";
    mkdir "$outputpath/2_Syn_NonSyn_filter";
    mkdir "$outputpath/3_Syn_NonSyn_filter_aa";
    open (CODINGRG,"$options{'codingRegion'}");
    open (CODINGRGR,">$outputpath/CodingRegion_r.txt");
    while (<CODINGRG>) {
        if (/^Protein	Beg/) {
            print CODINGRGR;
        }else{
            chomp;
            my @eachprebed=split(/\t/);
            print CODINGRGR "$eachprebed[0]\t$eachprebed[1]\t$eachprebed[2]\tconsensus\n";
        }
    }
    close CODINGRG; close CODINGRGR;

    system ("perl $selfpath/Scripts/diversiutils_modified_v2.pl -bam $outputpath/alignment/$list\_markup.bam -ref $outputpath/alignment/$list\.support_oem/consensus.txt -orfs $outputpath/CodingRegion_r.txt -stub $outputpath/1_Syn_NonSyn_aa/$list");
    system ("perl $selfpath/Scripts/diversifilter_v2.pl -in $outputpath/1_Syn_NonSyn_aa/$list -pQ 0.05 -pS 100000 -stub $outputpath/2_Syn_NonSyn_filter/$list");
    system ("perl $selfpath/Scripts/Syn_NonSyn_filter_aa_new.pl $outputpath/1_Syn_NonSyn_aa/$list\_entropy.txt $outputpath/CodingRegion_r.txt $outputpath/2_Syn_NonSyn_filter/$list\_filter.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_condon.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_AA.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_condon_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_AA_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_minor_change_filtered.txt");
    unlink ("$outputpath/CodingRegion_r.txt");
    &parseminorchange("$options{'codingRegion'}", "$outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_condon_filtered.txt","$outputpath/3_Syn_NonSyn_filter_aa/$list\_minor_change_filtered_ref.txt");
}

sub illuminaamplicon {
    ###mapping run 1
    print "The sample name is $list\n";
    my $indexfileExist = -e "$options{'ref'}\.1.bt2";
    if ($indexfileExist) {
        print "The bowtie2-build has been done\n";
    } else {
        print "The bowtie2-build is running\n";
        system ("bowtie2-build $options{'ref'} $options{'ref'}");
    }
    
    mkdir "$outputpath/alignment";
    if (exists $options{'fq'} and !exists $options{'fq1'} and !exists $options{'fq2'}) {
        if (exists $options{'maxins'}) {
            print "\"-maxins\" is only for the paired-end illumina ampilcon sequencing\n";
            exit;
        }
        print $options{'fq'}, "\n";
        system("bowtie2 --local -x $options{'ref'} -p $thread -U $options{'fq'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2308 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }else {
        if (!exists $options{'maxins'}) {
            print "Please provide maximum fragment length of in your amplicon Library\n";
            exit;
        }
        system("bowtie2 --local -X $options{'maxins'} --no-mixed -x $options{'ref'} -p $thread -1 $options{'fq1'} -2 $options{'fq2'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2316 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }
    
    system ("samtools index $outputpath/alignment/$list\.sorted.bam");
    #unlink ("$outputpath/alignment/$list\.sam");
    
    ###bam clipper run 1
    system ("bash $selfpath/Scripts/bamclipper.sh -b $outputpath/alignment/$list\.sorted.bam -p $outputpath/bamclipper_primer_bed.txt -n $thread");
    unlink ("$outputpath/alignment/$list\.sorted.bam");
    unlink ("$outputpath/alignment/$list\.sorted.bam.bai");
    system ("mv ./$list\.sorted.primerclipped.bam* $outputpath/alignment");
    system ("samtools view -h $outputpath/alignment/$list\.sorted.primerclipped.bam > $outputpath/alignment/$list\.sorted.primerclipped.sam");
    open(DATA,"$outputpath/alignment/$list\.sorted.primerclipped.sam");
    open(DATAR,">$outputpath/alignment/$list\.sorted.primerclipped.filtered.sam");
    my @readdata=<DATA>;
    close DATA;
    foreach (@readdata) {
        my @each=split(/\t/);
        if (/^@/) {
            print DATAR;
        }else {
            unless ($each[-4]=~/NM\:i/ and $each[-3]=~/MD\:Z/) {
                print DATAR;
            }  
        }
    }
    close DATA; close DATAR;
    system ("samtools view -@ $thread -Sb $outputpath/alignment/$list\.sorted.primerclipped.filtered.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam");
    system ("samtools index $outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam");
    system("java -XX:+UseParallelGC -XX:+UseNUMA -XX:NewRatio=9 -Xms10G -Xmx200G -jar $selfpath/Scripts/QuasiRecomb.jar -i $outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam -o $outputpath/alignment/$list\.support_oem -coverage");
    system("perl $selfpath/Scripts/parse_quasirecomb_V2_new.pl $outputpath/alignment/$list\.support_oem/support/allel_distribution_phred_weighted.txt $outputpath/alignment/$list\.support_oem/support/coverage.txt $outputpath/alignment/$list\.support_oem/consensus.txt");
    system ("rm -r $outputpath/alignment/$list\.support_oem/support");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.bam");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.bam.bai");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.sam");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.filtered.sam");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam.bai");
    
    ###mapping run 2
    system("bowtie2-build --threads $thread $outputpath/alignment/$list\.support_oem/consensus.txt $outputpath/alignment/$list\.support_oem/consensus");
    if (exists $options{'fq'} and !exists $options{'fq1'} and !exists $options{'fq2'}) {
        if (exists $options{'maxins'}) {
            print "\"-maxins\" is only for the paired-end illumina ampilcon sequencing\n";
            exit;
        }
        print $options{'fq'}, "\n";
        system("bowtie2 --local -x $outputpath/alignment/$list\.support_oem/consensus -p $thread -U $options{'fq'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2308 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }else {
        if (!exists $options{'maxins'}) {
            print "Please provide maximum fragment length of in your amplicon Library\n";
            exit;
        }
        system("bowtie2 --local -X $options{'maxins'} --no-mixed -x $outputpath/alignment/$list\.support_oem/consensus -p $thread -1 $options{'fq1'} -2 $options{'fq2'} -S $outputpath/alignment/$list\.sam 2> $outputpath/alignment/$list\_mapping_rate");
        system ("samtools view -@ $thread -q 10 -F 2316 -Sb $outputpath/alignment/$list\.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.bam");
    }
    system ("samtools index $outputpath/alignment/$list\.sorted.bam");
    
    ###bam clipper run 2
    open (PREBED,"$outputpath/bamclipper_primer_bed.txt");
    open (PREBEDR,">$outputpath/bamclipper_primer_bed_r.txt");
    while (<PREBED>) {
        my @eachprebed=split(/\t/,$_,2);
        print PREBEDR "consensus\t$eachprebed[1]";
    }
    close PREBED; close PREBEDR;

    system ("bash $selfpath/Scripts/bamclipper.sh -b $outputpath/alignment/$list\.sorted.bam -p $outputpath/bamclipper_primer_bed_r.txt -n $thread");
    unlink ("$outputpath/alignment/$list\.sorted.bam");
    unlink ("$outputpath/alignment/$list\.sorted.bam.bai");
    system ("mv ./$list\.sorted.primerclipped.bam* $outputpath/alignment");
    system ("samtools view -h $outputpath/alignment/$list\.sorted.primerclipped.bam > $outputpath/alignment/$list\.sorted.primerclipped.sam");
    open(DATA,"$outputpath/alignment/$list\.sorted.primerclipped.sam");
    open(DATAR,">$outputpath/alignment/$list\.sorted.primerclipped.filtered.sam");
    my @readdata=<DATA>;
    close DATA;
    foreach (@readdata) {
        my @each=split(/\t/);
        if (/^@/) {
            print DATAR;
        }else {
            unless ($each[-4]=~/NM\:i/ and $each[-3]=~/MD\:Z/) {
                print DATAR;
            }  
        }
    }
    close DATA; close DATAR;
    system ("samtools view -@ $thread -Sb $outputpath/alignment/$list\.sorted.primerclipped.filtered.sam | samtools sort -@ $thread -o $outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam");
    system ("samtools index $outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.bam");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.bam.bai");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.sam");
    unlink ("$outputpath/alignment/$list\.sorted.primerclipped.filtered.sam");

    ###running diversiutils, filtration, minor variation
    mkdir "$outputpath/1_Syn_NonSyn_aa";
    mkdir "$outputpath/2_Syn_NonSyn_filter";
    mkdir "$outputpath/3_Syn_NonSyn_filter_aa";
    open (CODINGRG,"$options{'codingRegion'}");
    open (CODINGRGR,">$outputpath/CodingRegion_r.txt");
    while (<CODINGRG>) {
        if (/^Protein	Beg/) {
            print CODINGRGR;
        }else{
            chomp;
            my @eachprebed=split(/\t/);
            print CODINGRGR "$eachprebed[0]\t$eachprebed[1]\t$eachprebed[2]\tconsensus\n";
        }
    }
    close CODINGRG; close CODINGRGR;

    system ("perl $selfpath/Scripts/diversiutils_modified_v2.pl -bam $outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam -ref $outputpath/alignment/$list\.support_oem/consensus.txt -orfs $outputpath/CodingRegion_r.txt -stub $outputpath/1_Syn_NonSyn_aa/$list");
    system ("perl $selfpath/Scripts/diversifilter_v2.pl -in $outputpath/1_Syn_NonSyn_aa/$list -pQ 0.05 -pS 100000 -stub $outputpath/2_Syn_NonSyn_filter/$list");
    system ("perl $selfpath/Scripts/Syn_NonSyn_filter_aa_new.pl $outputpath/1_Syn_NonSyn_aa/$list\_entropy.txt $outputpath/CodingRegion_r.txt $outputpath/2_Syn_NonSyn_filter/$list\_filter.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_condon.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_AA.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_condon_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_AA_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_minor_change_filtered.txt");
    unlink ("$outputpath/bamclipper_primer_bed_r.txt");
    unlink ("$outputpath/CodingRegion_r.txt");
    &parseminorchange("$options{'codingRegion'}", "$outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_condon_filtered.txt","$outputpath/3_Syn_NonSyn_filter_aa/$list\_minor_change_filtered_ref.txt");
}

sub nanoporeamplicon {
    if (!exists $options{'maxins'}) {
        print "Please provide maximum fragment length of in your amplicon Library\n";
        exit;
    }
    
    ###preparing fastq
    mkdir "$outputpath/alignment";
    mkdir "$outputpath/tmpfasta";
    if ($options{'fq'}=~/\.gz/) {
        system ("cp $options{'fq'} $outputpath/tmpfasta/tmpfasta.fastq.gz");
    }else{
        system ("cp $options{'fq'} $outputpath/tmpfasta/tmpfasta.fastq");
    }
    
    #system ("chmod 777 $selfpath/Scripts/artic");
    #system ("chmod 777 $selfpath/Scripts/align_trim");
    system ("artic guppyplex --skip-quality-check --min-length $minins --max-length $options{'maxins'} --directory $outputpath/tmpfasta --output $outputpath/alignment/lengthed_$list\.fastq");
    unlink ("$outputpath/tmpfasta");
    
    ###mapping
    system ("minimap2 -a -x map-ont -t $thread $options{'ref'} $outputpath/alignment/lengthed_$list\.fastq | samtools view -bS -F 4 - | samtools sort -o $outputpath/alignment/$list\.sorted.bam -");
    system ("align_trim --normalise 0 $options{'primerbed'} --remove-incorrect-pairs --report $outputpath/alignment/$list\.alignreport.txt < $outputpath/alignment/$list\.sorted.bam 2> $outputpath/alignment/$list\.alignreport.er | samtools sort -T $list - -o $outputpath/alignment/$list\.primertrimmed.rg.sorted.bam");
    system ("samtools index $outputpath/alignment/$list\.primertrimmed.rg.sorted.bam");
    unlink ("$outputpath/alignment/$list\.sorted.bam");
    
    ###running diversiutils, filtration, minor variation
    mkdir "$outputpath/1_Syn_NonSyn_aa";
    system ("perl $selfpath/Scripts/diversiutils_modified_v2.pl -bam $outputpath/alignment/$list\.primertrimmed.rg.sorted.bam -ref $options{'ref'} -orfs $options{'codingRegion'} -stub $outputpath/1_Syn_NonSyn_aa/$list");
}

sub primerbed {
    open(BED, "$options{'primerbed'}");
    my @bed=<BED>;
    chomp @bed;
    close BED;
    
    open(BEDR, ">$outputpath/bamclipper_primer_bed.txt");
    for (my $nbed=1; $nbed<$#bed+1; $nbed++) {
        my @leftprimers=grep (/\_$nbed\_LEFT/, @bed);
        my @rightprimers=grep (/\_$nbed\_RIGHT/, @bed);
    
        if ($#leftprimers == -1) {
            last;
        }
    
        for (my $n=0; $n<$#leftprimers+1; $n++) {
            my @eachlineleft=split(/\t/,$leftprimers[$n]);
            my @eachlineright=split(/\t/,$rightprimers[$n]);
            print BEDR "$eachlineleft[0]\t$eachlineleft[1]\t$eachlineleft[2]\t$eachlineright[0]\t$eachlineright[1]\t$eachlineright[2]\n";
        }
    }
    close BEDR;
}

sub parseminorchange {
    my %c2p;
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

    my @genome;
    open (GENOMESEQ,"$options{'ref'}");
    while (<GENOMESEQ>) {
        chomp;
        unless (/^\>/) {
        push (@genome,"$_");
        }
    }
    close GENOMESEQ;
    my $genomesequence=join("", @genome);
    my @eachgenomesequence=split(//,$genomesequence);

    open (CODINGREG, "$options{'codingRegion'}");
    open (SVGAA, ">$outputpath/soruce_virus_genome_aa.txt");
    print SVGAA "Protein\tAAPosition\tRefAA\tRefSite\tRefCodon\n";
    while (<CODINGREG>) {
        unless (/^Protein\tBeg/) {
            my @begend=split(/\t/);
            my @orig=@eachgenomesequence[$begend[1]-1 .. $begend[2]-1];
            my @vars;
            push @vars, [ splice @orig, 0, 3 ] while @orig;

            my $AAPosition=0;
            my $refSite=$begend[1]-3;
            foreach my $var(@vars) {
                $AAPosition++;
                $refSite=$refSite+3;
                print SVGAA "$begend[0]\t$AAPosition\t",$c2p{uc(join("",@{$var}))},"\t$refSite\t",@{$var}, "\n";
           }
        }
    }
    close CODINGREG; close SVGAA;

    my ($newinputREGION, $newinputCONDONALLR,$newoutputCONDONALLR)=($_[0],$_[1],$_[2]);
    print "$newinputREGION, $newinputCONDONALLR,$newoutputCONDONALLR";
    open (GENOME,"$outputpath/soruce_virus_genome_aa.txt");
    my %hashgenoment; my %hashgenomeaa;
    while (<GENOME>) {
        chomp;
        my @eachgenome=split(/\t/);
        $hashgenoment{"$eachgenome[0]\t$eachgenome[1]"}=$eachgenome[4];
        $hashgenomeaa{"$eachgenome[0]\t$eachgenome[1]"}=$eachgenome[2];
    }
    close GENOME;

    my %hashelecments;
    open (CONDONALLR, "$newinputCONDONALLR");
    while (<CONDONALLR>){
        chomp;
        unless (/^Chr	Protein/) {
            my @eachelecments=split(/\t/);
            my @groupelecments=("$eachelecments[1]","$eachelecments[2]","$eachelecments[3]",$hashgenomeaa{"$eachelecments[1]\t$eachelecments[2]"},"$eachelecments[5]",$hashgenoment{"$eachelecments[1]\t$eachelecments[2]"},"$eachelecments[7]","$eachelecments[8]","$eachelecments[9]");
            my $newaa=$hashgenomeaa{"$eachelecments[1]\t$eachelecments[2]"};
            push (@{$hashelecments{"$eachelecments[1]\t$eachelecments[2]\t$newaa\t$eachelecments[9]"}}, [@groupelecments]);
        }
    }
    close CONDONALLR;

    open(REGION, "$newinputREGION");
    my @regionlist=<REGION>;
    shift @regionlist;
    close REGION;

    open(MINORCHANGE, ">$newoutputCONDONALLR");
    print MINORCHANGE "Protein\tAAPosition\tRefAA\tAAcoverage\tCntNonSyn\tCntSyn\tCntStop\tCntN\n";
    foreach (@regionlist) {
        chomp;
        my @eachcodongregions=split(/\t/);

        foreach my $key1 (sort {$hashelecments{$a}[0][1] <=> $hashelecments{$b}[0][1]} keys %hashelecments) {
            if ($key1=~/^$eachcodongregions[0]\t/) {
                #print "$key1 => ";
                print MINORCHANGE "$key1\t";
                my @aaref=split(/\t/, $key1);
                my @codoncollection;
                my %hashaacount;
                foreach my $elementingroup(@{$hashelecments{"$key1"}}) {
                    unless (@$elementingroup[5] eq @$elementingroup[6]) {
                        #print "@$elementingroup[5]\:@$elementingroup[6]";
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
                    #print "\:\:\: nonsym\:", sum(@nonsym),"\ ";
                    print MINORCHANGE sum(@nonsym),"\t";
                }else{
                    #print "\:\:\: nonsym\:", 0 ,"\ ";
                    print MINORCHANGE "0\t";
                }
    
                if ($#sym >= 0) {
                    #print " sym\:", sum(@sym),"\ ";
                    print MINORCHANGE sum(@sym),"\t";
                }else{
                    #print " sym\:", 0 ,"\ ";
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
    close CONDONALLR; close MINORCHANGE;
}


