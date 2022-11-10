# LycSpaceTimePoolSeq
Notes on my analyses of our pool-seq Lycaeides space-time data set. See [Pool_DNA_plates](https://drive.google.com/drive/folders/1U4AsshyMvlySNtODuSLWo0dYDH_rgho4) for details on samples and sample concentrations.

# Data
The pool-seq data are currently in /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/. The current list of samples follows (this is from a single round of sequencing):

| Population | Samples (N)|
|------------|------------|
| ABM | ABM20 (48) |
| BCR | BCR17 (48), BCR17rep (48) |
| BHP | BHP19 (48) |
| BKM | BKM19 (33) |
| BTB | BTB10 (38), BTB13 (48), BTB13rep (48), BTB17 (48), BTB17rep (48) |
| CLH | CLH19 (36) |
| CP | CP98 (37), CP01 (44), CP19 (48), CP19rep (48) |
| EP | EP19 (48), EP19rep (48) | 
| GNP | GNP17 (56) |
| HJ | HJ20 (48) |
| HNV | HNV17 (48) |
| LS | LS19 (48) |
| MEN | MEN12 (10) |
| MR | MR98 (44), MR20 (48) |
| MTU | MTU20 (48) |
| SBW | SBW18 (20) |
| SHC | SHC11 (46) |
| SIN | SIN10 (48) |
| SUV | SUV20 (51) |
| TBY | TBY11 (24) |
| TIC | TIC19 (48) |
| VE | VE06 (48), VE20 (48) |
| YG | YT98 (33), YG20 (48) |

Replicates (rep) are the same DNA extractions re-pooled and sequenced as a distinct sample. MEN = outgroup. The data were generated and cleaned up by BGI. Reports for each round of sequencing are attached: [BGI_F22FTSUSAT0310-01_LYCgpswR_report_en.pdf](https://github.com/zgompert/LycSpaceTimePoolSeq/files/9940314/BGI_F22FTSUSAT0310-01_LYCgpswR_report_en.pdf)
. 

# DNA Sequence Alignment

I am aligning the DNA sequence data to the updated (based on PacBio) *L. melissa* genome. I am using `bwa-mem2` for this, which is basically just a sped up version of `bwa mem` that also works directly with gzipped files [https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2). I am using `bwa-mem2` version 2.0pre2. 

First, I (re)indexed the reference genome.

```bash
## index genome with bwa-mem2
/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta
```
Then, I set up the alignment. The submission scrip is (from /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts):

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bwa-mem2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignments

perl BwaMemFork.pl ../F22FTSUSAT0310-01_LYCgpswR/soapnuke/clean/*/*1.fq.gz 
```

Which runs the following:

```perl
#!/usr/bin/perl
#
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";

FILES:
foreach $fq1 (@ARGV){
	$pm->start and next FILES; ## fork
	$fq2 = $fq1;
	$fq2 =~ s/_1\.fq\.gz/_2.fq.gz/ or die "failed substitution for $fq1\n";
        $fq1 =~ m/clean\/([A-Za-z0-9]+)/ or die "failed to match id $fq1\n";
	$ind = $1;
	$fq1 =~ m/([A-Za-z_\-0-9]+)_1\.fq\.gz$/ or die "failed match for file $fq1\n";
	$file = $1;
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:Lyc-"."$ind\\tLB:Lyc-"."$ind\\tSM:Lyc-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$file.bam - && samtools index -@ 2 $ind"."_$file.bam\n";

	$pm->finish;
}
```
I am pipping the results on to `samtools` (version 1.16) to compress, sort and index the alignments.

I then merged the bam files for each population using `samtools` version 1.16. This was donw with the following shell script:

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=merge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignment

perl ../Scripts/MergeFork.pl 
```

Which runs


```perl
#!/usr/bin/perl
#
# merge alignments for each population sample with samtools version XX 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

open(IDS,"pids.txt");
while(<IDS>){
	chomp;
	push(@IDs,$_);
}
close(IDS);

FILES:
foreach $id (@IDs){
	$pm->start and next FILES; ## fork
        system "samtools merge -c -p -o Merged/$id.bam $id"."_*.bam\n";
	system "samtools index -@ 2 Merged/$id.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```

pids.txt lists all of the population IDs.
