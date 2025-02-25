</p>
AioMinor was implemented in Perl programming language, including a main script for generation of consensus (dominant) virus genome, and calling of nucletide and amino acid varation frequency in both dominant and minor level. It accepts bascalled fastq files derived from Nanopore amplicon sequencing, cleaned/trimmed Illumina amplicon/normal fastq files (single-end or paired-end) with a SARS-CoV-2/other virus genome. By default, AioMinor analyses SARS-CoV-2 based on an NCBI reference genome (NC_045512.2), but the user can also provide customized SARS-CoV-2 or other virus genome as a reference.

## Installation:
**1. Create an environment with one step**
```
git clone https://github.com/Hiscox-lab/AioMinor/AioMinor.git
cd AioMinor
conda env create -f my_environment.yml
source activate AioMinor 
```
**2. Create an environment step by step**

Third party dependencies:
  > samtools(>=1.9)
  > 
  > bowtie2(=2.4.1)
  > 
  > minimap2(=2.24)

All these third party tool dependencies should be exported to PATH, so that AioMinor can find them. 

Perl module dependencies:
```
Getopt::Long
File::Basename
Data::Dumper
IO::File
Math::CDF
List::Util
Text::NSP
Bio::DB::Sam
Bio::SeqIO
```

## Usage

Please see the details of each parameter by:

```
Required options:
  -platform           sequencing platform "nanopore" or "illumina".
  -method             method of sequencing Library Preparation "amplicon" or "cDNA".
  -maxins             maximum fragment length of in your amplicon Library with paired-end sequencing.
  -ref                reference genome sequence in fasta file.
  -codingRegion       coding regions in your reference genome.
  -primerbed          primer scheme in bed file.
  -fq                 fastq(.gz) file for single-end sequencing.
  -fq1                fastq(.gz) file for paired-end sequencing mate 1s.
  -fq2                fastq(.gz) file for paired-end sequencing mate 2s.
  -primerbed          amplicon primer information in bed file.

Optional options:
  -samplename         sample name, "alignment" by default.
  -maxins             maximum fragment length of in your amplicon sequencing Library,
  -minins             minimum fragment length of in your amplicon sequencing library, "50" by default.
  -rotation           number of the consensus genome polish, "2" by default.
  -t/-thread          number of threads, 1 by default.
  -o                  output path, "./" by default.
  
  -v/-version         Print version information.
  -h/-help            Print help message.\n\n";
```

## **Examples:**
**1. ARTIC Illumina sequencing data**

To analyse ARTIC-Illumina amplicon sequencing data with your chosen primer scheme. (The ARTIC primer bed can be found in the "Primerbeds" folder, and reference genome and coding region of SARS-CoV-2 NC_045512.2 can be found in the "References" folder):

```
perl AioMinor.pl -t 16 -platform illumina -method amplicon -maxins 500 -ref genome.fasta -codingRegion CodingRegion.txt -primerbed nCoV-2019.primer_V3.bed -fq1 example_R1.fastq.gz -fq2 example_R2.fastq.gz -o example_output
```

**2. ARTIC Nanopore sequencing data**

To analyse ARTIC-Nanopore amplicon sequencing data with your chosen primer scheme. (The ARTIC primer bed can be found in the "Primerbeds" folder, and reference genome and coding region of SARS-CoV-2 NC_045512.2 can be found in the "References" folder):
```
perl AioMinor.pl -t 16 -platform nanopore -method amplicon -maxins 1600 -ref genome.fasta -codingRegion CodingRegion.txt -primerbed nCoV-2019.primer_V3.bed -fq example.fastq.gz -o example_output
```

**3. Normal Illumina sequencing data**

To analyse Normal Illumina sequencing data with your chosen primer scheme. (The reference genome and coding region of SARS-CoV-2 NC_045512.2 can be found in the "References" folder):
```
perl AioMinor.pl -t 16 -platform nanopore -method amplicon -maxins 1600 -ref genome.fasta -codingRegion CodingRegion.txt -fq1 example_R1.fastq.gz -fq2 example_R2.fastq.gz -o example_output
```

## **Results**
The results can be found under the in output path.

### 1_Syn_NonSyn_aa**

*_entropy.txt file in this folder prvides the raw nucletide varation frequency compared to the consensus genome,inlcuding Transitions and transversions. *_AA.txt file in this folder prvides the raw amino acid varation frequency compared to the consensus genome,inlcuding synonymous and non-synonymous substitution. 

### 2_Syn_NonSyn_filter**

*_filter.txt file in this folder prvides the filtered nucletide varation frequency,inlcuding Transitions and transversions.

### 3_Syn_NonSyn_filter_aa**

*_AA_all_AA_filtered.txt file in this folder prvides the details of minor and major amino acids at each amino acid site after filtration in 2_Syn_NonSyn_filter. *_AA_all_condon_filtered.txt file in this folder prvides the details of minor and major condons at each amino acid site after filtration in 2_Syn_NonSyn_filter. *_minor_change_filtered.txt file in this folder prvides minor frequency of synonymous and non-synonymous substitution compared to the consensus genome after filtration. 

### alignment/alignment.support_oem**
consensus.txt file in this folder prvides consensus genome sequences.


## **Test data**
There a test data obtained from ARTIC (V3) Illumina sequencing of a cell culture sample in the Testdata folder. AioMinor can be tested with this data in the AioMinor directory.

```
perl AioMinor.pl -t 16 -platform illumina -method amplicon -maxins 500 -ref ./References/genome_NC_045512.2.fasta -codingRegion ./References/CodingRegion_NC_045512.2.txt -primerbed ./Primerbeds/nCoV-2019.primer_V3.bed -fq1 ./Testdata/cell_illumina_R1.fastq.gz -fq2 ./Testdata/cell_illumina_R2.fastq.gz -o cell_illumina_output
```

## Citations
Dong, Xiaofeng, et al. "Using minor variant genomes and machine learning to study the genome biology of SARS-CoV-2 over time." Nucleic Acids Research 53.4 (2025): gkaf077.

Dong, Xiaofeng, et al. "Variation around the dominant viral genome sequence contributes to viral load and outcome in patients with Ebola virus disease." Genome biology 21 (2020): 1-20.




