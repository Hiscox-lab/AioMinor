#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

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
    ###mapping
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
    
    ###running diversiutils, filtration, minor variation
    mkdir "$outputpath/1_Syn_NonSyn_aa";
    mkdir "$outputpath/2_Syn_NonSyn_filter";
    mkdir "$outputpath/3_Syn_NonSyn_filter_aa";
    system ("perl $selfpath/Scripts/diversiutils_modified_v2.pl -bam $outputpath/alignment/$list\_markup.bam -ref $options{'ref'} -orfs $options{'codingRegion'} -stub $outputpath/1_Syn_NonSyn_aa/$list");
    system ("perl $selfpath/Scripts/diversifilter_v2.pl -in $outputpath/1_Syn_NonSyn_aa/$list -pQ 0.05 -pS 100000 -stub $outputpath/2_Syn_NonSyn_filter/$list");
    system ("perl $selfpath/Scripts/Syn_NonSyn_filter_aa_new.pl $outputpath/1_Syn_NonSyn_aa/$list\_entropy.txt $options{'codingRegion'} $outputpath/2_Syn_NonSyn_filter/$list\_filter.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_condon.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_AA.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_condon_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_AA_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_minor_change_filtered.txt");
}

sub illuminaamplicon {
    ###mapping
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
    
    ###bam clipper
    system ("bash $selfpath/Scripts/bamclipper.sh -b $outputpath/alignment/$list\.sorted.bam -p $outputpath/bamclipper_primer_bed.txt -n $thread");
    unlink ("$outputpath/alignment/$list\.sorted.bam");
    #unlink ("$outputpath/alignment/$list\.sorted.bam.bai");
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
    #unlink ("$outputpath/alignment/$list\.sorted.primerclipped.bam");
    #unlink ("$outputpath/alignment/$list\.sorted.primerclipped.bam.bai");
    #unlink ("$outputpath/alignment/$list\.sorted.primerclipped.sam");
    #unlink ("$outputpath/alignment/$list\.sorted.primerclipped.filtered.sam");
    
    ###running diversiutils, filtration, minor variation
    mkdir "$outputpath/1_Syn_NonSyn_aa";
    mkdir "$outputpath/2_Syn_NonSyn_filter";
    mkdir "$outputpath/3_Syn_NonSyn_filter_aa";
    system ("perl $selfpath/Scripts/diversiutils_modified_v2.pl -bam $outputpath/alignment/$list\.sorted.primerclipped.filtered.sorted.bam -ref $options{'ref'} -orfs $options{'codingRegion'} -stub $outputpath/1_Syn_NonSyn_aa/$list");
    system ("perl $selfpath/Scripts/diversifilter_v2.pl -in $outputpath/1_Syn_NonSyn_aa/$list -pQ 0.05 -pS 100000 -stub $outputpath/2_Syn_NonSyn_filter/$list");
    system ("perl $selfpath/Scripts/Syn_NonSyn_filter_aa_new.pl $outputpath/1_Syn_NonSyn_aa/$list\_entropy.txt $options{'codingRegion'} $outputpath/2_Syn_NonSyn_filter/$list\_filter.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_condon.txt $outputpath/1_Syn_NonSyn_aa/$list\_AA_all_AA.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_condon_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_AA_all_AA_filtered.txt $outputpath/3_Syn_NonSyn_filter_aa/$list\_minor_change_filtered.txt");
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
    
    system ("chmod 777 $selfpath/Scripts/artic");
    system ("chmod 777 $selfpath/Scripts/align_trim");
    system ("$selfpath/Scripts/artic guppyplex --skip-quality-check --min-length $minins --max-length $options{'maxins'} --directory $outputpath/tmpfasta --output $outputpath/alignment/lengthed_$list\.fastq");
    unlink ("$outputpath/tmpfasta");
    
    ###mapping
    system ("minimap2 -a -x map-ont -t $thread $options{'ref'} $outputpath/alignment/lengthed_$list\.fastq | samtools view -bS -F 4 - | samtools sort -o $outputpath/alignment/$list\.sorted.bam -");
    system ("$selfpath/Scripts/align_trim --normalise 0 $options{'primerbed'} --remove-incorrect-pairs --report $outputpath/alignment/$list\.alignreport.txt < $outputpath/alignment/$list\.sorted.bam 2> $outputpath/alignment/$list\.alignreport.er | samtools sort -T $list - -o $outputpath/alignment/$list\.primertrimmed.rg.sorted.bam");
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


