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

I am aligning the DNA sequence data to the updated (based on PacBio) *L. melissa* genome. I am using `bwa-mem2` for this, which is basically just a sped up version of `bwa mem` that also works directly with gzipped files.

