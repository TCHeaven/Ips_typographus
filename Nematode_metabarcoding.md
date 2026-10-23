# Reanalysis of V. Rau <i>Ips typographus</i> / nematode 18S and 28S metabarcoding dataset

Veronica Rau originally generated a dataset of 18S and 28S metabarcoding reads as part of a project to investigate, primarily, the diversity of nematodes associated with <i>Ipy typographus</i> in the S. Tyrol region using molecular methods, as previous identification has been done via morphology only.

D2A and D3B primers were used to target the 28S region: ACAAGTACCGTGAGGGAAAGTTG + TCGGAAGGAACCAGCTACTA
NEM (forward and reverse) primers were used to target the 18S region: GCAAGTCTGGTGCCAGCAGC + CCGTGTTGAGTCAAATTAAG

These primers were chosen as they target the nematode 18 and 28S regions. However, excluding reads aligning to the host <i>I. typographus</i> genome, >25% of reads were classified to nematoda via BLASTP search, and <50% in any given sample. This left many reads unassigned.

We are interested to know what other species are represented in the dataset amongst the non-nematoda, non-Ips reads. Maja Fluch investigated further by using the BOLD (Barcode Of Life Data Systems) database to assign identity to the reads as an alternative approach to BLASTP against the NCBI nt database (Veronica's method). However, the results of this analysis remain confusing as many taxa are assigned that seem biologically implausible.

My aim is the draw out some clearer patterns from the data - if possible. Below are documented analysis and commands used with this dataset. All commands were executed on the ScientificNet HPC Cluster accessed at the University of Bozen-Bolzano from the directory /data/users/theaven/nematode_project, unless stated otherwise.

## Contents

1. [Collecting Data](#2)
2. [Quality Control and ASV Inference](#3)
  2.1 [FastQC](#7)
  2.2 [Cutadapt](#8)
  2.3 [Fastp](#9)
  2.4 [Bowtie2](#10)
  2.5 [DADA2](#11)
  2.6 [Minimum-abundance / prevalence filtering](#12)
3. [Taxonomy Assignment - QIIME2](#6)
  3.1 [QIIME2](#17)
  3.2 [QIIME2 Custom](#18)
4. [Taxonomy Assignment - IDTAXA](#4)
  4.1 [QIIME2](#15)
  4.2 [IDTAXA](#16)
  4.3 [IDTAXA Custom](#19)
5. [Taxonomy Assignment - BLASTN](#5)
  5.1 [BLASTN](#14)
  5.2 [NCBI nt](#20)
  5.3 [SILVA SSU and LSU - BLAST locally / ACT lookup](#21)
6. [Taxonomy Assignment - BOLD](#22)
7. [Nematode Curated Databases](#23)
  7.1 [IDTAXA - nemabiome](#24)
  7.2 [IDTAXA - nemabase](#25)
8. [Relate my ASVs to Veronika/Maja's reads](#26)

## Collecting data <a name="2"></a>
Raw sequencing reads were retreived from the archive folder \\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nematode_metabarcoding and uploaded to the HPC:
```bash
ls /data/users/theaven/nematode_project/raw_data/18S
ls /data/users/theaven/nematode_project/raw_data/28S

for file in $(ls /data/users/theaven/nematode_project/raw_data/*/*.fastq.gz); do
ID=$(basename $file | rev | cut -d '_' -f2- | rev)
mkdir $(dirname $file)/$ID
mv $file $(dirname $file)/$ID/.
done
```

## Quality Control and ASV Inference <a name="3"></a>

### FastQC  <a name="7"></a>
The raw sequence reads were subjected to a quality control check using FastQC.
```bash
for ReadDir in /data/users/theaven/nematode_project/raw_data/*S/*; do
    ID=$(basename "$ReadDir")
    echo "$ID"
    zcat "$ReadDir"/*1.fastq.gz | wc -l | awk '{print $1/4}'
    zcat "$ReadDir"/*2.fastq.gz | wc -l | awk '{print $1/4}'
done

screen -S nematode
module load anaconda3
for ReadDir in $(ls -d /data/users/theaven/nematode_project/raw_data/*S/*); do
	Task=FastQC
	ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
	OutDir="$ReadDir"/"$Task"
	ExpectedOutput="$OutDir"/$(basename "${Reads[0]}" | sed 's@.fastq.gz@@g')_fastqc.html

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 5s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_fastqc.sh "$OutDir" "${Reads[@]}")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done
```
### Cutadapt  <a name="8"></a>
Primers were removed from the reads where present using Cutadapt.

NOTE:The reads are a mix of paired and single end samples.
```bash
screen -r nematode
for ReadDir in $(ls -d /data/users/theaven/nematode_project/raw_data/28S/*); do
	Task=CutAdapt
	ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
	OutDir="$(echo "$ReadDir" | sed 's@raw_data@qc_data@g')/"$Task""
	Forward_Primer=ACAAGTACCGTGAGGGAAAGTTG
	Reverse_Primer=TCGGAAGGAACCAGCTACTA
	ExpectedOutput="$OutDir"/$(basename "${Reads[0]}" | sed 's@.fastq.gz@.trim.fastq.gz@g')

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 5s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_cutadapt.sh "$OutDir" "$Forward_Primer" "$Reverse_Primer" "${Reads[@]}")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done

for ReadDir in $(ls -d /data/users/theaven/nematode_project/raw_data/18S/*); do
	Task=CutAdapt
	ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
	OutDir="$(echo "$ReadDir" | sed 's@raw_data@qc_data@g')/"$Task""
	echo "$OutDir"
	Forward_Primer=GCAAGTCTGGTGCCAGCAGC
	Reverse_Primer=CCGTGTTGAGTCAAATTAAG
	ExpectedOutput="$OutDir"/$(basename "${Reads[0]}" | sed 's@.fastq.gz@.trim.fastq.gz@g')

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 5s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_cutadapt.sh "$OutDir" "$Forward_Primer" "$Reverse_Primer" "${Reads[@]}")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done

for ReadDir in /data/users/theaven/nematode_project/qc_data/*/*/CutAdapt; do
    ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    echo "$ID" >> post_cutadapt_count.txt
    zcat "$ReadDir"/*1.trim.fastq.gz | wc -l | awk '{print $1/4}' >> post_cutadapt_count.txt
    zcat "$ReadDir"/*2.trim.fastq.gz | wc -l | awk '{print $1/4}' >> post_cutadapt_count.txt
done
```
### Fastp  <a name="9"></a>

Reads were trimmed with Fastp to remove adapters, reads/pairs shorter than 100bp or with >40% of bases below phred 20 were discarded.
```bash
screen -r nematode
for ReadDir in $(ls -d /data/users/theaven/nematode_project/qc_data/*S/*/CutAdapt); do
	Task=Fastp
	ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
	OutDir="$(dirname "$ReadDir")/"$Task""
	ExpectedOutput="$OutDir"/$(basename "${Reads[0]}" | sed 's@.fastq.gz@.trimmed.fastq.gz@g')

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 60s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_fastp.sh "$OutDir" "${Reads[@]}")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done

for ReadDir in /data/users/theaven/nematode_project/qc_data/*/*/Fastp; do
    ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    echo "$ID" >> post_fastp_count.txt
    zcat "$ReadDir"/*1.trim.trimmed.fastq.gz | wc -l | awk '{print $1/4}' >> post_fastp_count.txt
    zcat "$ReadDir"/*2.trim.trimmed.fastq.gz | wc -l | awk '{print $1/4}' >> post_fastp_count.txt
done

#For 28S samples 10s of thousands of reads are retained after fastp, for 18S there are some samples where read count is reduced below 10,000 (still >3,000). I will continue with fastp filtered reads as this step was performed by V. Rau.
```
### Bowtie2  <a name="10"></a>

Reads were aligned to the reference <i>I. typographus</i> genome in order to exclude reads from the host.
```bash
screen -r nematode
srun -p cpu  -c 4 --mem 16G --pty bash
module load anaconda3 
conda activate bowtie2
cd ~/genomes/Ips/typographus/GCA_016097725.1
bowtie2-build --threads 4 GCA_016097725.1_CZU_Ityp_1.0_genomic.fna GCA_016097725.1_CZU_Ityp_1.0_genomic_index
exit()

for ReadDir in $(ls -d /data/users/theaven/nematode_project/qc_data/*S/*/Fastp); do
	Task=Bowtie2
	ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
	OutDir="$(dirname "$ReadDir")/"$Task""
	Reference=~/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic_index
	ExpectedOutput="${OutDir}"/"${ID}"_vs_"$(basename $Reference | sed 's@_index@@g')"_mapped.bam

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 300s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_bowtie2.sh "$OutDir" "$ID" "$Reference" "${Reads[@]}")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done

#Remove intermediates:
rm -r /data/users/theaven/nematode_project/qc_data/*S/*/CutAdapt
rm /data/users/theaven/nematode_project/qc_data/*S/*/Fastp/*.fastq.gz
```
### DADA2  <a name="11"></a>
The denoiser tool DADA2 was run to model and remove error patterns from the Illumina data: read ends where quality drops were trimmed, reads with high expected error (EE), or above a max EE threshold, or with any ambiguous bases (N) were discarded. DADA2’s removeBimeraDenovo function was also used to remove chimera from the ASV table. For paired reads DADA2 also merges forward and reverse reads in the overlapping region. When F and R disagree, the merge algorithm uses quality scores to pick the most likely base or discards the read. The output of DADA2 is an Alternative Sequence Variant (ASV) table as well as filtered denoised read fastq files.

DADA2 is an R package, 18S and 28S, paired and single end reads must be run seperately. The relavent files were downloaded:
```bash
for Dir in $(ls -d /data/users/theaven/nematode_project/qc_data/*S/*/Bowtie2); do
    if ls "$Dir"/*.1.fastq.gz 1> /dev/null 2>&1 && ls "$Dir"/*.2.fastq.gz 1> /dev/null 2>&1; then
        Out=/data/users/theaven/download_18122025/$(echo "$Dir" | cut -d '/' -f7)/paired/$(echo "$Dir" | cut -d '/' -f8)
        mkdir -p "$Out"
        cp "$Dir"/*unaligned*.fastq.gz "$Out"/.
    else
        Out=/data/users/theaven/download_18122025/$(echo "$Dir" | cut -d '/' -f7)/single/$(echo "$Dir" | cut -d '/' -f8)
        mkdir -p "$Out"
        cp "$Dir"/*unaligned*.fastq.gz "$Out"/.
    fi
done
```
"/data/users/theaven/download_18122025" was subseqeuntly downloaded to "C:\Users\THeaven\OneDrive - Scientific Network South Tyrol\R"

**Plot reads and select appropriate truncation lengths:**

DADA2 filter settings will truncate reads, a parameter is required below which length the read is discarded entirely. Plots of read quality help to determine a reasonable cuttoff, ie. before major dropoff in quality of many reads, or the median read quality is < Q25-30, or the lower percentile < Q20, mean and median should also be similar. However, shorter, more permissive lengths are prefered as other filter settings should act as a safety net to remove problem reads. However, for paired reads the truncated forward and reverse reads must still be long enough to overlap (at least 20bp) given the amplicon size or the reads cannot be merged.
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("dada2","ShortRead", "Biostrings"))
install.packages("dplyr")

library(dada2)
library(ShortRead)
library(Biostrings)
library(dplyr)

setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")
set.seed(1)

#18S
#Single end example reads plotted:
plotQualityProfile("download_18122025/18S/single/NEM_ST2_8/18S_NEM_ST2_8_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #30-180 
plotQualityProfile("download_18122025/18S/single/NEM_ST2_9/18S_NEM_ST2_9_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #110-200-250 
plotQualityProfile("download_18122025/18S/single/NEM_RV061/18S_NEM_RV061_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #140-225-275 
plotQualityProfile("download_18122025/18S/single/NEM_RV062/18S_NEM_RV062_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #140-225-275 
plotQualityProfile("download_18122025/18S/single/NEM_RV106/18S_NEM_RV106_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #175-260-275 

#28S
#Single end example reads plotted:
plotQualityProfile("download_18122025/28S/single/D2A_RV061/28S_D2A_RV061_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #175-220-250 
plotQualityProfile("download_18122025/28S/single/D2A_RV062/28S_D2A_RV062_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #160-220-250 
plotQualityProfile("download_18122025/28S/single/D2A_ST2_8/28S_D2A_ST2_8_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #160-160-200 
plotQualityProfile("download_18122025/28S/single/D2A_ST2_9/28S_D2A_ST2_9_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #175-200-225 
plotQualityProfile("download_18122025/28S/single/D2A_RV082/28S_D2A_RV082_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.fastq.gz") #175-200-250 

#Paired end example reads plotted:
plotQualityProfile("download_18122025/28S/paired/RV224/28S_RV224_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.1.fastq.gz") #200-250-275 
plotQualityProfile("download_18122025/28S/paired/RV224/28S_RV224_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.2.fastq.gz") #140-190-210 

plotQualityProfile("download_18122025/28S/paired/RV129/28S_RV129_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.1.fastq.gz") #median starts to decline ~200-260-275 
plotQualityProfile("download_18122025/28S/paired/RV129/28S_RV129_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.2.fastq.gz") #140-190-200 

plotQualityProfile("download_18122025/28S/paired/RU1_1/28S_RU1_1_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.1.fastq.gz") #220-260-275 
plotQualityProfile("download_18122025/28S/paired/RU1_1/28S_RU1_1_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.2.fastq.gz") #150-90-220 

plotQualityProfile("download_18122025/28S/paired/CR1_9/28S_CR1_9_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.1.fastq.gz") #190-260-275 
plotQualityProfile("download_18122025/28S/paired/CR1_9/28S_CR1_9_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.2.fastq.gz") #150-190-210 

plotQualityProfile("download_18122025/28S/paired/Ab1_10/28S_Ab1_10_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.1.fastq.gz") #200-260-275 
plotQualityProfile("download_18122025/28S/paired/Ab1_10/28S_Ab1_10_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.2.fastq.gz") #140-190-210 
```
DADA2 quality profile plots summarise per-cycle base quality scores across reads. The solid orange line is median quality score at each base position, the solid turquoise line is mean quality score at each position, orange dashed lines show 10th and 90th percentiles. Q20 = ~1% error rate. Q30 = ~ 0.1% error rate.

![Representative DADA2 quality profile plot](figures/28S_D2A_RV082_unaligned_vs_GCA_016097725.1_CZU_Ityp_1.0_genomic.png)

For subset of samples plotted the mean (turquoise) is frequently noticably lower than the median (orange), this indicates that the quality distribution is skewed: many reads are decent (keeping the median high), but a subset of reads are really bad (dragging the mean down). maxEE is the primary setting to target these problem reads.

**Collect inputs to run DADA2:**
```R
single_18S <- list.files(path = "download_18122025/18S/single", full.names = TRUE, recursive = TRUE)
single_28S <- list.files(path = "download_18122025/28S/single", full.names = TRUE, recursive = TRUE)
paired_28S <- list.files(path = "download_18122025/28S/paired", full.names = TRUE, recursive = TRUE)

single_18S_1 <- single_18S[!grepl("\\.[12]\\.fastq\\.gz$", single_18S)]
single_28S_1 <- single_28S[!grepl("\\.[12]\\.fastq\\.gz$", single_28S)]
paired_28S_1 <- paired_28S[grepl("\\.1\\.fastq.gz$", paired_28S)]
paired_28S_2 <- paired_28S[grepl("\\.2\\.fastq.gz$", paired_28S)]

get_samplename <- function(x) sub("_unaligned_.*$", "", basename(x))

single_18S_lookup_table <- data.frame(
  sample = sapply(single_18S_1, get_samplename),
  fn = single_18S_1,
  stringsAsFactors = FALSE
)

single_28S_lookup_table <- data.frame(
  sample = sapply(single_28S_1, get_samplename),
  fn = single_28S_1,
  stringsAsFactors = FALSE
)

snF <- vapply(paired_28S_1, get_samplename, character(1))
snR <- vapply(paired_28S_2, get_samplename, character(1))
paired_samples <- intersect(snF, snR)
paired_28S_lookup_table <- data.frame(
  sample = paired_samples,
  fnF = paired_28S_1[match(paired_samples, snF)],
  fnR = paired_28S_2[match(paired_samples, snR)],
  stringsAsFactors = FALSE
)
```
**Run DADA2 for the different datasets:**
truncLen was selected based upon where median median base quality score starts to drop in the plotted samples. A stricter maxEE setting was selected to try to filter out problem reads of particularly poor quality.

Run DADA2 for 18S single end reads:
```R
#Parameters
truncLen <- 110  #Truncate reads to this length, below this remove entirely
truncQ <- 20     #Truncate at first base with quality <= truncQ (higher = more conservative | lower = more permissive)
maxEE <- 1       #Maximum expected errors (1 is strict, 2 is moderate stringency, higher = more sensitivity + more errors tolerated)
maxN <- 0        #Maximum allowed N bases
rm.phix <- TRUE  #Remove PhiX reads (bacteriophage used as control in Illumina sequencing runs)
pool <- FALSE    #Pool samples for error rate estimation (FALSE = error model is learned per sample, this is more conservative and faster/less memory | pseudo = max sensitivity)
threads <- TRUE  #Use multithreading
outdir <- "download_18122025/results/single_18S_dada2_output"
in_df <- single_18S_lookup_table

#Create output directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Filter and trim reads
filt_paths <- file.path(outdir, paste0(in_df$sample, "_filt.fastq.gz"))
out <- filterAndTrim(fwd = in_df$fn, filt = filt_paths, truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, rm.phix = rm.phix, compress = TRUE, multithread = threads)

#Learn error rates, dereplicate, and run DADA2 algorithm
err <- learnErrors(filt_paths, nbases = 5e10, multithread = threads)
derep <- derepFastq(filt_paths, verbose = TRUE)
names(derep) <- in_df$sample
dada_out <- dada(derep, err = err, pool = pool, multithread = threads)

#Make ASV table and remove chimeras
seqtab <- makeSequenceTable(dada_out)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = TRUE)
nonchim_reads <- rowSums(seqtab.nochim)
chim_reads    <- rowSums(seqtab) - nonchim_reads
chim_track <- cbind(nonchim = nonchim_reads, chimera = chim_reads, chimera_fraction = chim_reads / rowSums(seqtab))


#Save results
saveRDS(seqtab.nochim, file = file.path(outdir, "seqtab.nochim.rds"))
write.table(seqtab.nochim, file.path(outdir, "seqtab.nochim.tsv"), sep="\t", quote=FALSE, col.names=NA)
dim(seqtab.nochim)
#  80 1247

getN <- function(x) sum(getUniques(x))
track <- cbind(input = out[, "reads.in"], filtered = out[, "reads.out"], denoised = sapply(dada_out, getN))
rownames(track) <- in_df$sample
write.table(track, file.path(outdir, "read_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)

write.table(chim_track, file.path(outdir, "chimera_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
```
52.4% of single end 18S reads were retained through filtering and denoising. Of these filtered reads, 0 - 6.7% were subseqeuently identified as chimeric depending on the sample, the median portion of chimeric reads was 1%. 1247 ASVs were inferred.

Run DADA2 for 28S single end reads:
```R
#Parameters
truncLen <- 160  #Truncate reads to this length, then remove entirely
truncQ <- 20     #Truncate at first base with quality <= truncQ (higher = more conservative | lower = more permissive)
maxEE <- 1       #Maximum expected errors (1 is strict, 2 is moderate stringency, higher = more sensitivity + more errors tolerated)
maxN <- 0        #Maximum allowed N bases
rm.phix <- TRUE  #Remove PhiX reads (bacteriophage used as control in Illumina sequencing runs)
pool <- FALSE    #Pool samples for error rate estimation (FALSE = error model is learned per sample, this is more conservative and faster/less memory | pseudo = max sensitivity)
threads <- TRUE  #Use multithreading
outdir <- "download_18122025/results/single_28S_dada2_output"
in_df <- single_28S_lookup_table

#Create output directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Filter and trim reads
filt_paths <- file.path(outdir, paste0(in_df$sample, "_filt.fastq.gz"))
out <- filterAndTrim(fwd = in_df$fn, filt = filt_paths, truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, rm.phix = rm.phix, compress = TRUE, multithread = threads)

#Learn error rates, dereplicate, and run DADA2 algorithm
err <- learnErrors(filt_paths, nbases = 5e10, multithread = threads)
derep <- derepFastq(filt_paths, verbose = TRUE)
names(derep) <- in_df$sample
dada_out <- dada(derep, err = err, pool = pool, multithread = threads)

#Make ASV table and remove chimeras
seqtab <- makeSequenceTable(dada_out)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = TRUE)
nonchim_reads <- rowSums(seqtab.nochim)
chim_reads    <- rowSums(seqtab) - nonchim_reads
chim_track <- cbind(nonchim = nonchim_reads, chimera = chim_reads, chimera_fraction = chim_reads / rowSums(seqtab))


#Save results
saveRDS(seqtab.nochim, file = file.path(outdir, "seqtab.nochim.rds"))
write.table(seqtab.nochim, file.path(outdir, "seqtab.nochim.tsv"), sep="\t", quote=FALSE, col.names=NA)

getN <- function(x) sum(getUniques(x))
track <- cbind(input = out[, "reads.in"], filtered = out[, "reads.out"], denoised = sapply(dada_out, getN))
rownames(track) <- in_df$sample
write.table(track, file.path(outdir, "read_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)

write.table(chim_track, file.path(outdir, "chimera_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
```
54.5% of single end 28S reads were retained through filtering and denoising. Of these filtered reads, 0 - 11.6% were subseqeuently identified as chimeric depending on the sample, the median portion of chimeric reads was 3.8%.

Run DADA2 for paired end reads:
```R
#Parameters
truncLen <- c(190, 140)  #Truncate reads to this length, then remove entirely
truncQ <- 20     #Truncate at first base with quality <= truncQ (higher = more conservative | lower = more permissive)
maxEE <- c(1,3)  #Maximum expected errors (4 for the reverse reads as these are known/expected to be lower quality)
maxN <- 0        #Maximum allowed N bases
rm.phix <- TRUE  #Remove PhiX reads (bacteriophage used as control in Illumina sequencing runs)
pool <- FALSE    #Pool samples for error rate estimation (FALSE = error model is learned per sample, this is more conservative and faster/less memory | pseudo = max sensitivity)
threads <- TRUE  #Use multithreading
outdir <- "download_18122025/results/paired_28S_dada2_output"
in_df <- paired_28S_lookup_table

#Create output directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Filter and trim reads (paired)
filtFs <- file.path(outdir, paste0(in_df$sample, "_F_filt.fastq.gz"))
filtRs <- file.path(outdir, paste0(in_df$sample, "_R_filt.fastq.gz"))
out <- filterAndTrim(fwd = in_df$fnF, filt = filtFs, rev = in_df$fnR, filt.rev = filtRs, truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, rm.phix = rm.phix, compress = TRUE, multithread = threads)

#Learn error rates, dereplicate, run DADA2 algorithm, and merge pairs
#Learn errors separately for F and R
errF <- learnErrors(filtFs, nbases = 5e10, multithread = threads)
errR <- learnErrors(filtRs, nbases = 5e10, multithread = threads)
#Derep separately
derepF <- derepFastq(filtFs, verbose = TRUE)
derepR <- derepFastq(filtRs, verbose = TRUE)
names(derepF) <- in_df$sample
names(derepR) <- in_df$sample
#Denoise separately
dadaF <- dada(derepF, err = errF, pool = pool, multithread = threads)
dadaR <- dada(derepR, err = errR, pool = pool, multithread = threads)
#Merge pairs
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)

#Make ASV table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = TRUE)
nonchim_reads <- rowSums(seqtab.nochim)
chim_reads    <- rowSums(seqtab) - nonchim_reads
chim_track <- cbind(nonchim = nonchim_reads, chimera = chim_reads, chimera_fraction = chim_reads / rowSums(seqtab))

# Save results
saveRDS(seqtab.nochim, file = file.path(outdir, "seqtab.nochim.rds"))
write.table(seqtab.nochim, file.path(outdir, "seqtab.nochim.tsv"), sep="\t", quote=FALSE, col.names=NA)

getN <- function(x) sum(getUniques(x))
track <- cbind(input    = out[, "reads.in"], filtered = out[, "reads.out"], denoisedF = sapply(dadaF, getN), denoisedR = sapply(dadaR, getN), merged = sapply(mergers, function(x) sum(x$abundance)))
rownames(track) <- in_df$sample
write.table(track, file.path(outdir, "read_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)

write.table(chim_track, file.path(outdir, "chimera_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
```
Very few read pairs could be merged successfully, even when very permissive filtering settings were trialled (truncQ = 1, maxEE = 2,4, truncLen = 100,100), this is likely because the reads generated from Illumina sequencing are 300bp in length and the 28S amplicon is expected to be 580-650bp depending on the species, this leaves no margin for trimming of primers and lower quality bases, indeed in most cases even with no trimming at all paired reads will not have the 20bp minimum overlap required for merging.

The higher quality forward reads from the 28S paired data were combined with the single end 28S reads to produce a combined D2A forward primer dataset. It is clear from the previous analyses that the quality of the single end reads is lower than the paired forward reads, therefore different error models are calculated for each; however, the same truncLen is used for both to harmonise the ASVs - this means a more permissive (160) setting for the paired forward reads:
```R
#Parameters
truncLen <- 160  #Truncate reads to this length, or remove entirely if below this length
truncQ <- 20     #Truncate at first base with quality <= truncQ (higher = more conservative | lower = more permissive)
maxEE <- 1       #Maximum expected errors (1 is strict, 2 is moderate stringency, higher = more sensitivity + more errors tolerated)
maxN <- 0        #Maximum allowed N bases
rm.phix <- TRUE  #Remove PhiX reads (bacteriophage used as control in Illumina sequencing runs)
pool <- FALSE    #Pool samples for error rate estimation (FALSE = error model is learned per sample, this is more conservative and faster/less memory | pseudo = max sensitivity)
threads <- TRUE  #Use multithreading
outdir <- "download_18122025/results/combined_28S_dada2_output"
in_df <- paired_28S_lookup_table

#Create output directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Filter and trim reads
filt_paths <- file.path(outdir, paste0(in_df$sample, "_filt.fastq.gz"))
out <- filterAndTrim(fwd = in_df$fnF, filt = filt_paths, truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, rm.phix = rm.phix, compress = TRUE, multithread = threads)

#Learn error rates, dereplicate, and run DADA2 algorithm
err <- learnErrors(filt_paths, nbases = 5e10, multithread = threads)
derep <- derepFastq(filt_paths, verbose = TRUE)
names(derep) <- in_df$sample
dada_out <- dada(derep, err = err, pool = pool, multithread = threads)

#Make ASV table and remove chimeras
seqtab <- makeSequenceTable(dada_out)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = TRUE)
nonchim_reads <- rowSums(seqtab.nochim)
chim_reads    <- rowSums(seqtab) - nonchim_reads
chim_track <- cbind(nonchim = nonchim_reads, chimera = chim_reads, chimera_fraction = chim_reads / rowSums(seqtab))

#Combine
seqtab_single_nochim <- readRDS("download_18122025/results/single_28S_dada2_output/seqtab.nochim.rds")
seqtab_28S_all_nochim <- mergeSequenceTables(seqtab_single_nochim, seqtab.nochim)
dim(seqtab_28S_all_nochim)
#  254 3322
dim(seqtab_single_nochim)
#   80 1785
dim(seqtab.nochim)
#  174 1649

getN <- function(x) sum(getUniques(x))
track <- cbind(input = out[, "reads.in"], filtered = out[, "reads.out"], denoised = sapply(dada_out, getN))
rownames(track) <- in_df$sample
track_single <- read.table("download_18122025/results/single_28S_dada2_output/read_tracking_truncQ-20_maxEE-1_truncLen-160.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
common_cols <- intersect(colnames(track_single), colnames(track))
track_all <- rbind(track_single[, common_cols, drop = FALSE], track[, common_cols, drop = FALSE])

chim_single <- read.table("download_18122025/results/single_28S_dada2_output/chimera_tracking.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
common_cols_chim <- intersect(colnames(chim_single), colnames(chim_track))
chim_all <- rbind(chim_single[, common_cols_chim, drop = FALSE], chim_track[,  common_cols_chim, drop = FALSE])

#Save results
saveRDS(seqtab.nochim, file = file.path(outdir, "seqtab.nochim.rds"))
write.table(seqtab.nochim, file.path(outdir, "seqtab.nochim.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(track, file.path(outdir, "read_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(chim_track, file.path(outdir, "chimera_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)

saveRDS(seqtab_28S_all_nochim, file = file.path(outdir, "combined_seqtab.nochim.rds"))
write.table(seqtab_28S_all_nochim, file.path(outdir, "combined_seqtab.nochim.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(track_all, file.path(outdir, "combined_read_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(chim_all, file.path(outdir, "combined_chimera_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
```
68.4% of paired forward 28S reads were retained through filtering and denoising. Of these filtered reads, 0 - 9.7% were subseqeuently identified as chimeric depending on the sample, the median portion of chimeric reads was 1%. ASVs resulting from chimeric reads were removed from the ASV table (these and denoise reads remain in the filtered reads fastq.gz files). 1785 ASVs were inferred from the single end 28S dataset, 1649 ASVs were inferred from the paired end dataset. When combining the single end 28S and paired end 28S forward reads only 100 of 3322 ASVs overlapped between the two groups - although the overlap in taxa represented in these ASVs may be greater...

### Minimum-abundance / prevalence filtering  <a name="12"></a>

ASVs that occur in only one sample and represent less than 0.1% of reads from that sample were dropped. (Sequencing errors and tag-jumps almost always show up in one sample only)
```R
ASV_18S <- readRDS("download_18122025/results/single_18S_dada2_output/seqtab.nochim.rds")
ASV_28S <- readRDS("download_18122025/results/combined_28S_dada2_output/combined_seqtab.nochim.rds")

prune_single_sample_low_abundance <- function(seqtab, min_rel = 0.001) {
  prevalence <- colSums(seqtab > 0) #count how many samples each ASV appears in
  sample_totals <- rowSums(seqtab) #find total reads per sample
  rel_abund <- sweep(seqtab, 1, sample_totals, "/") #ASV_reads / total_reads_in_sample
  rel_abund[is.na(rel_abund)] <- 0 #Set relative abundance to 0 where undefined
  max_rel <- apply(rel_abund, 2, max) #Finds the highest relative abundance ASV ever reaches (for single sample ASVs there is only 1 value)
  keep <- !(prevalence == 1 & max_rel < min_rel) #remove ASV if only in one sample and highest abundance is below setting
  seqtab[, keep, drop = FALSE]
}

ASV_18S_pruned <- prune_single_sample_low_abundance(ASV_18S, min_rel = 0.001)
ASV_28S_pruned <- prune_single_sample_low_abundance(ASV_28S, min_rel = 0.001)

ncol(ASV_18S) #1247
ncol(ASV_18S_pruned) #629
mean(rowSums(ASV_18S_pruned > 0)) #33.6

ncol(ASV_28S) #3322
ncol(ASV_28S_pruned) #1604
mean(rowSums(ASV_28S_pruned > 0)) #20.0

dir.create("download_18122025/ASV_18S", showWarnings = FALSE, recursive = TRUE)
dir.create("download_18122025/ASV_28S", showWarnings = FALSE, recursive = TRUE)

saveRDS(ASV_18S_pruned, file = file.path("download_18122025/ASV_18S/ASV_18S_asv.rds"))
write.table(ASV_18S_pruned, file.path("download_18122025/ASV_18S/ASV_18S_asv.tsv"), sep="\t", quote=FALSE, col.names=NA)
saveRDS(ASV_28S_pruned, file = file.path("download_18122025/ASV_28S/ASV_28S_asv.rds"))
write.table(ASV_28S_pruned, file.path("download_18122025/ASV_28S/ASV_28S_asv.tsv"), sep="\t", quote=FALSE, col.names=NA)

# seqtab_filt: samples x ASVs
asv_seqs <- colnames(ASV_18S_pruned)
asv_headers <- paste0("ASV", seq_along(asv_seqs))
dna <- Biostrings::DNAStringSet(asv_seqs)
names(dna) <- asv_headers
Biostrings::writeXStringSet(dna, "download_18122025/ASV_18S/ASVs.fasta")
write.csv(data.frame(ASV=asv_headers, Sequence=asv_seqs), "download_18122025/ASV_18S/ASV_id_map.csv", row.names = FALSE)

asv_seqs <- colnames(ASV_28S_pruned)
asv_headers <- paste0("ASV", seq_along(asv_seqs))
dna <- Biostrings::DNAStringSet(asv_seqs)
names(dna) <- asv_headers
Biostrings::writeXStringSet(dna, "download_18122025/ASV_28S/ASVs.fasta")
write.csv(data.frame(ASV=asv_headers, Sequence=asv_seqs), "download_18122025/ASV_28S/ASV_id_map.csv", row.names = FALSE)
````
## Taxonomy Assignment - QIIME2 <a name="6"></a>

For high-throughput classification the most commonly used appraoch in the literature appears to be QIIME2 classifier + SILVA database. QIIME2 is a naive Bayes classifier that learns k-mer frequency patterns, it handles short, conserved sequences better than BLAST and is better for assigning ambiguous reads, handling subtle variations in conserved regions, and avoiding overconfident species matches. The output is probability-based taxonomic assignments. 

Prepare inputs DADA2(R) -> QIIME2(HPC):
```bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

ASV_18S_dir=/data/users/theaven/nematode_project/asvs/ASV_18S
ASV_28S_dir=/data/users/theaven/nematode_project/asvs/ASV_28S

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/make_qiime2_inputs.py \
  --counts "$ASV_18S_dir"/ASV_18S_asv.tsv \
  --map "$ASV_18S_dir"/ASV_id_map.csv \
  --out-dir "$ASV_18S_dir"/qiime_inputs \
  --out-prefix qiime_inputs

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/make_qiime2_inputs.py \
  --counts "$ASV_28S_dir"/ASV_28S_asv.tsv \
  --map "$ASV_28S_dir"/ASV_id_map.csv \
  --out-dir "$ASV_28S_dir"/qiime_inputs \
  --out-prefix qiime_inputs
```

### QIIME2 <a name="17"></a>

Pre-trained QIIME2 classifiers for SSU and LSU regions are available from SILVA. However, these are prepared on the full length region, not the amplicon produced by our specific set of primers. For high-level taxonomy (domain, phylum, maybe class) this should be fine to use. For genus- or species-level resolution accuracy will likely be low.

Visualise data with QIIME:
```bash
apptainer exec docker://quay.io/qiime2/amplicon:2025.7 qiime --version
apptainer pull ~/git_repos/Containers/qiime2-amplicon-2025.7.sif docker://quay.io/qiime2/amplicon:2025.7

#Visualise inputs:
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i "$ASV_18S_dir"/qiime_inputs/qiime_inputsASV_table.tsv \
  -o "$ASV_18S_dir"/qiime_inputs/ASV_table.biom \
  --table-type="OTU table" \
  --to-hdf5

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i "$ASV_28S_dir"/qiime_inputs/qiime_inputsASV_table.tsv \
  -o "$ASV_28S_dir"/qiime_inputs/ASV_table.biom \
  --table-type="OTU table" \
  --to-hdf5

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path "$ASV_18S_dir"/qiime_inputs/ASV_table.biom \
  --output-path "$ASV_18S_dir"/qiime_inputs/table.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path "$ASV_28S_dir"/qiime_inputs/ASV_table.biom \
  --output-path "$ASV_28S_dir"/qiime_inputs/table.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table summarize \
  --i-table "$ASV_18S_dir"/qiime_inputs/table.qza \
  --o-visualization "$ASV_18S_dir"/qiime_inputs/table.qzv

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table summarize \
  --i-table "$ASV_28S_dir"/qiime_inputs/table.qza \
  --o-visualization "$ASV_28S_dir"/qiime_inputs/table.qzv

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table tabulate-seqs \
  --i-data "$ASV_18S_dir"/rep-seqs.qza \
  --o-visualization "$ASV_18S_dir"/rep-seqs.qzv

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table tabulate-seqs \
  --i-data "$ASV_28S_dir"/rep-seqs.qza \
  --o-visualization "$ASV_28S_dir"/rep-seqs.qzv
```

Run QIIME with pre-trained full lenth SSU and LSU classifiers:
```bash
srun -p cpu  -c 8 --mem 32G --pty bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

ASV_18S_dir=/data/users/theaven/nematode_project/asvs/ASV_18S
ASV_28S_dir=/data/users/theaven/nematode_project/asvs/ASV_28S

#Import ASVs
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$ASV_18S_dir"/ASVs.fasta \
  --output-path "$ASV_18S_dir"/rep-seqs.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$ASV_28S_dir"/ASVs.fasta \
  --output-path "$ASV_28S_dir"/rep-seqs.qza

#Classify
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-classifier classify-sklearn \
  --i-classifier ~/db/SILVA/ssu_138.2/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza \
  --i-reads "$ASV_18S_dir"/rep-seqs.qza \
  --o-classification "$ASV_18S_dir"/taxonomy.qza \
  --p-n-jobs 8

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-classifier classify-sklearn \
  --i-classifier ~/db/SILVA/lsu_138.2/SILVA138.2_LSURef_NR99_uniform_classifier_full-length.qza \
  --i-reads "$ASV_28S_dir"/rep-seqs.qza \
  --o-classification "$ASV_28S_dir"/taxonomy.qza \
  --p-n-jobs 8

#Export
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_18S_dir"/taxonomy.qza \
  --output-path "$ASV_18S_dir"/qiime_fullength_taxonomy

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_28S_dir"/taxonomy.qza \
  --output-path "$ASV_28S_dir"/qiime_fullength_taxonomy

#Visualise:

#I do not have the metadata for the samples (location, year, etc.), therefore create a minimal metadata file:
echo -e "#SampleID" > "$ASV_18S_dir"/sample-metadata.tsv
head -n 1 "$ASV_18S_dir"/qiime_inputs/qiime_inputsASV_table.tsv | cut -f2- | tr '\t' '\n' >> "$ASV_18S_dir"/sample-metadata.tsv

echo -e "#SampleID" > "$ASV_28S_dir"/sample-metadata.tsv
head -n 1 "$ASV_28S_dir"/qiime_inputs/qiime_inputsASV_table.tsv | cut -f2- | tr '\t' '\n' >> "$ASV_28S_dir"/sample-metadata.tsv

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa barplot \
  --i-table "$ASV_18S_dir"/qiime_inputs/table.qza \
  --i-taxonomy "$ASV_18S_dir"/taxonomy.qza \
  --m-metadata-file "$ASV_18S_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_18S_dir"/taxa-barplot.qzv

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa barplot \
  --i-table "$ASV_28S_dir"/qiime_inputs/table.qza \
  --i-taxonomy "$ASV_28S_dir"/taxonomy.qza \
  --m-metadata-file "$ASV_28S_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_28S_dir"/taxa-barplot.qzv
```

### QIIME2 Custom <a name="18"></a>
Whilst pre-trained QIIME2 classifiers for SSU and LSU regions are available from SILVA these are prepared on the full length region, or for different primers to ours. QIIME2 Naïve Bayes ideally needs a database trimmed to your exact primer region. For ASVs < ~500 bp (typical amplicons) region-specific training is usually worth it. Our ASVs are 110 and 160bp - much shorter than the full region. It is therefore advised to train a classifier on SILVA trimmed to our primer region to hopefully get a per-ASV taxonomy up to genus/species (if present in the DB). 

A new classifier was trained for the NEM-F 18S ASVs. The NEM forward primer is used with a very permissive reverse primer and a truncation length equal to the filtered reads/ASVs:
```bash
#Import SILVA SSU reference files into QIIME2
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[RNASequence]' \
  --input-path ~/db/SILVA/ssu_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta \
  --output-path ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-seqs.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[SILVATaxidMap]' \
  --input-path ~/db/SILVA/ssu_138.2/taxmap_slv_ssu_ref_nr_138.2.txt \
  --output-path ~/db/SILVA/ssu_138.2/taxmap-silva-138.2-ssu-nr99.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[SILVATaxonomy]' \
  --input-path ~/db/SILVA/ssu_138.2/tax_slv_ssu_138.2.txt \
  --output-path ~/db/SILVA/ssu_138.2/taxranks-silva-138.2-ssu-nr99.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'Phylogeny[Rooted]' \
  --input-path ~/db/SILVA/ssu_138.2/tax_slv_ssu_138.2.tre \
  --output-path ~/db/SILVA/ssu_138.2/taxtree-silva-138.2-nr99.qza

#Generate a fixed‑rank taxonomy with RESCRIPt
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime rescript parse-silva-taxonomy \
  --i-taxonomy-tree ~/db/SILVA/ssu_138.2/taxtree-silva-138.2-nr99.qza \
  --i-taxonomy-map ~/db/SILVA/ssu_138.2/taxmap-silva-138.2-ssu-nr99.qza \
  --i-taxonomy-ranks ~/db/SILVA/ssu_138.2/taxranks-silva-138.2-ssu-nr99.qza \
  --p-ranks domain kingdom phylum class order family genus \
  --p-include-species-labels \
  --p-no-rank-propagation \
  --o-taxonomy ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-taxonomy.qza

#Export
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-taxonomy.qza \
  --output-path ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-taxonomy

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-seqs.qza \
  --output-path ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-seqs

#Convert to DNA
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime rescript reverse-transcribe \
  --i-rna-sequences ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-seqs.qza \
  --o-dna-sequences ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-seqs-dna.qza

#Filter for amplicon region
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif \
qiime feature-classifier extract-reads \
  --i-sequences ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-seqs-dna.qza \
  --p-f-primer GCAAGTCTGGTGCCAGCAGC \
  --p-r-primer A \
  --p-read-orientation both \
  --p-trunc-len 110 \
  --o-reads ~/db/SILVA/ssu_138.2/NEM-F-reference-amplicon-110.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table tabulate-seqs \
  --i-data ~/db/SILVA/ssu_138.2/NEM-F-reference-amplicon-110.qza \
  --o-visualization ~/db/SILVA/ssu_138.2/NEM-F-reference-amplicon-110.qzv
#49441 seqeunces after primer filtering and truncation

#Train classifier
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ~/db/SILVA/ssu_138.2/NEM-F-reference-amplicon-110.qza \
  --i-reference-taxonomy ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-taxonomy.qza \
  --o-classifier ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-NEM-F-classifier.qza
```

A new classifier was trained for the D2A 28S ASVs. The D2A forward primer is used with a very permissive reverse primer and a truncation length equal to the filtered reads/ASVs. For the 28S region SILVA SSU reference files were already imported into QIIME2 as part of IDTAXA classification (4.1).
```bash
#Convert to DNA
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime rescript reverse-transcribe \
  --i-rna-sequences ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-seqs.qza \
  --o-dna-sequences ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-seqs-dna.qza

#Filter for amplicon region
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif \
qiime feature-classifier extract-reads \
  --i-sequences ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-seqs-dna.qza \
  --p-f-primer ACAAGTACCGTGAGGGAAAGTTG \
  --p-r-primer A \
  --p-read-orientation both \
  --p-trunc-len 160 \
  --o-reads ~/db/SILVA/lsu_138.2/D2A-F-reference-amplicon-160.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table tabulate-seqs \
  --i-data ~/db/SILVA/lsu_138.2/D2A-F-reference-amplicon-160.qza \
  --o-visualization ~/db/SILVA/lsu_138.2/D2A-F-reference-amplicon-160.qzv
#4662 seqeunces after primer filtering and truncation

#Train classifier
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ~/db/SILVA/lsu_138.2/D2A-F-reference-amplicon-160.qza \
  --i-reference-taxonomy ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-taxonomy.qza \
  --o-classifier ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-D2A-F-classifier.qza
```
Run QIIME2 classification with NEM-F and D2A-F custom classifiers
```bash
ASV_18S_dir=/data/users/theaven/nematode_project/asvs/ASV_18S
ASV_28S_dir=/data/users/theaven/nematode_project/asvs/ASV_28S

#Classify
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-classifier classify-sklearn \
  --i-classifier ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-NEM-F-classifier.qza \
  --i-reads "$ASV_18S_dir"/rep-seqs.qza \
  --o-classification "$ASV_18S_dir"/NEM-F-taxonomy.qza \
  --p-n-jobs 8

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-classifier classify-sklearn \
  --i-classifier ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-D2A-F-classifier.qza \
  --i-reads "$ASV_28S_dir"/rep-seqs.qza \
  --o-classification "$ASV_28S_dir"/D2A-F-taxonomy.qza \
  --p-n-jobs 8

#Export
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_18S_dir"/NEM-F-taxonomy.qza \
  --output-path "$ASV_18S_dir"/qiime_NEM-F_taxonomy

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_28S_dir"/D2A-F-taxonomy.qza \
  --output-path "$ASV_28S_dir"/qiime_D2A-F_taxonomy

#Visualise:
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa barplot \
  --i-table "$ASV_18S_dir"/qiime_inputs/table.qza \
  --i-taxonomy "$ASV_18S_dir"/NEM-F-taxonomy.qza \
  --m-metadata-file "$ASV_18S_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_18S_dir"/NEM-F-taxa-barplot.qzv

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa barplot \
  --i-table "$ASV_28S_dir"/qiime_inputs/table.qza \
  --i-taxonomy "$ASV_28S_dir"/D2A-F-taxonomy.qza \
  --m-metadata-file "$ASV_28S_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_28S_dir"/D2A-F-taxa-barplot.qzv
```
## Taxonomy Assignment - IDTAXA <a name="4"></a>

A ready-made IDTAXA training set is publically available for SILVA SSU and was downloaded from https://www2.decipher.codes/Downloads.html - Accessed 07/01/2026

No ready-made IDTAXA training set is publically available for SILVA LSU. A training set was built from the SILVA files downloaded from https://www.arb-silva.de/current-release/Exports - Accessed 07/01/2026

### QIIME2 <a name="16"></a>

In order to build an LSU training set the vaiable toxonomic information assigned to different entries in the SILVA database must first be converted into fixed-rank taxonomy information for each entry. QIIME 2’s parse-silva-taxonomy was used to do this:
```bash
module load qiime2/2024.10

#Import SILVA LSU reference files into QIIME2
qiime tools import \
  --type 'FeatureData[RNASequence]' \
  --input-path ~/db/SILVA/lsu_138.2/SILVA_138.2_LSURef_NR99_tax_silva.fasta \
  --output-path ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-seqs.qza

qiime tools import \
  --type 'FeatureData[SILVATaxidMap]' \
  --input-path ~/db/SILVA/lsu_138.2/taxmap_slv_lsu_ref_nr_138.2.txt \
  --output-path ~/db/SILVA/lsu_138.2/taxmap-silva-138.2-lsu-nr99.qza

qiime tools import \
  --type 'FeatureData[SILVATaxonomy]' \
  --input-path ~/db/SILVA/lsu_138.2/tax_slv_lsu_138.2.txt \
  --output-path ~/db/SILVA/lsu_138.2/taxranks-silva-138.2-lsu-nr99.qza

qiime tools import \
  --type 'Phylogeny[Rooted]' \
  --input-path ~/db/SILVA/lsu_138.2/tax_slv_lsu_138.2.tre \
  --output-path ~/db/SILVA/lsu_138.2/taxtree-silva-138.2-nr99.qza

#Generate a fixed‑rank taxonomy with RESCRIPt
qiime rescript parse-silva-taxonomy \
  --i-taxonomy-tree ~/db/SILVA/lsu_138.2/taxtree-silva-138.2-nr99.qza \
  --i-taxonomy-map ~/db/SILVA/lsu_138.2/taxmap-silva-138.2-lsu-nr99.qza \
  --i-taxonomy-ranks ~/db/SILVA/lsu_138.2/taxranks-silva-138.2-lsu-nr99.qza \
  --p-ranks domain kingdom phylum class order family genus \
  --p-include-species-labels \
  --p-no-rank-propagation \
  --o-taxonomy ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-taxonomy.qza

#Export
qiime tools export \
  --input-path ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-taxonomy.qza \
  --output-path ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-taxonomy

qiime tools export \
  --input-path ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-seqs.qza \
  --output-path ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-seqs
```
### IDTAXA <a name="13"></a>

IDTAXA (DECIPHER package in R) - probabilistic sequence classification, model-based learning - usually more accurate and conservative than BLAST or naive Bayes. IDTAXA does not require trimming and handles full-length sequences correctly. IDTAXA is slower than QIIME2 NB, but usually much more accurate for 18S/28S metabarcoding. If trained on full-length, IDTAXA often backs off to a safer rank rather than confidently guessing too deep.

Finished building the SILVA LSU IDTAXA training set using DECIPHER’s LearnTaxa() in R:

NOTE: the 'species' data in SILVA can be questionable as users list host organism in this field or write things like 'metagenome'. We will need to be carefull when interpreting the Assignment results at the species level. Also, for some entries species is given but not all the higher level taxa, this can cause problems because IDTAXA likes unique species to have only one taxonomic path. Other entries like Incertae Sedis also had to be made unambiguous. I have therefore changed every column to give the full ':' seperated taxonomic path to that level in order to bypass this problem, this tranformation will presumably need to be reversed in the outputs of IDTAXA in order to make them readable. 
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DECIPHER","stringr", "Biostrings"))

library(DECIPHER)
library(Biostrings)
library(stringr)

setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")

#Import data
tax_table <- read.delim("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/download_07012026/silva-138.2-lsu-nr99-taxonomy/taxonomy.tsv",
                        sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
ref_seqs_rna <- readRNAStringSet("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/download_07012026/silva-138.2-lsu-nr99-seqs/rna-sequences.fasta")
ref_seqs <- as(ref_seqs_rna, "DNAStringSet")
taxid <- read.delim("tax_slv_lsu_138.2.txt", header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE,
                    col.names = c("hierarchy", "Index", "Rank", "tmp", "version"))
seq_id <- sub(" .*", "", names(ref_seqs))
m <- match(seq_id, tax_table$Feature.ID)
tax_table2 <- tax_table[m, ]
stopifnot(all(tax_table2$Feature.ID == seq_id))
all(seq_id == tax_table2$Feature.ID) #should return 'TRUE'

#Clean and prepare taxonomy strings
tax_clean <- str_replace_all(tax_table2$Taxon, "[a-z]__","")
tax_list <- str_split(tax_clean, ";\\s*")
tax_mat <- do.call(rbind, tax_list)
colnames(tax_mat) <- c("domain","kingdom","phylum","class", "order","family","genus","species")
tax_mat[tax_mat == ""] <- NA

propagate_missing <- function(row) {
  last_known <- NA_character_
  for (i in seq_along(row)) {
    if (!is.na(row[i])) {
      last_known <- row[i]
    } else if (!is.na(last_known)) {
      row[i] <- paste0(last_known, "_unclassified")
    }
  }
  row
}

tax_mat <- t(apply(tax_mat, 1, propagate_missing))

expand_tax_mat_to_paths <- function(tax_mat, sep = ":") {
  tax_mat <- as.matrix(tax_mat)
  n_rows <- nrow(tax_mat)
  n_cols <- ncol(tax_mat)
  out <- tax_mat
  for (j in seq_len(n_cols)) {
    out[, j] <- apply(
      tax_mat[, seq_len(j), drop = FALSE],
      1,
      function(r) {
        r <- r[!is.na(r)]
        if (length(r) == 0) NA_character_
        else paste(r, collapse = sep)
      }
    )
  }
  colnames(out) <- colnames(tax_mat)
  out
}

#all levels
tax_mat2 <- expand_tax_mat_to_paths(tax_mat)
tax_strings2 <- apply(tax_mat2, 1, function(r) {paste(c("Root", r[!is.na(r)]), collapse="; ")})
tax_strings2 <- trimws(tax_strings2)
tax_strings2 <- gsub("\\s*;\\s*", "; ", tax_strings2)

#species only
tax_mat[, "species"] <- apply(tax_mat, 1, function(r) paste(r[!is.na(r)], collapse=":"))
tax_strings <- apply(tax_mat, 1, function(r) {paste(c("Root", r[!is.na(r)]), collapse="; ")})
tax_strings <- trimws(tax_strings)
tax_strings <- gsub("\\s*;\\s*", "; ", tax_strings)

#Build rank taxid table from the taxonomy strings
level_to_rank <- c("root", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
spl <- strsplit(tax_strings2, ";\\s*")
edge_list <- do.call(rbind, lapply(spl, function(x) {
  x <- x[x != ""]
  paths <- vapply(seq_along(x), function(i) paste(x[1:i], collapse="; "), character(1))
  parents <- c(NA_character_, paths[-length(paths)])
  data.frame(
    path = paths,
    parent_path = parents,
    level = seq_along(paths) - 1L,           # Root = 0
    name = x,
    stringsAsFactors = FALSE
  )
}))
nodes <- unique(edge_list[, c("path", "parent_path", "level", "name")])
nodes <- nodes[order(nodes$level, nodes$path), ]
nodes$Index <- seq_len(nrow(nodes)) - 1L
path_to_index <- setNames(nodes$Index, nodes$path)
nodes$Parent <- ifelse(is.na(nodes$parent_path), -1L, path_to_index[nodes$parent_path])
nodes$Rank <- level_to_rank[pmin(nodes$level + 1L, length(level_to_rank))]
rank_table3 <- nodes[, c("Index", "name", "Parent", "level", "Rank")]
colnames(rank_table3) <- c("Index", "Name", "Parent", "Level", "Rank")

#Check for any remaining ambiguous names
ambiguous <- rank_table3[
  duplicated(rank_table3[, c("Name", "Level")]) |
  duplicated(rank_table3[, c("Name", "Level")], fromLast = TRUE),
  c("Name", "Level")
]
unique(ambiguous)

#Train the IDTAXA classifier
#With rank
trainingSet3 <- LearnTaxa(
  train    = ref_seqs,
  taxonomy = tax_strings2,
  rank     = rank_table3,
  K        = NULL,       #NULL lets DECIPHER choose k‑mer size automatically
  maxIterations = 10,	 #Number of training iterations
  verbose  = TRUE
)
nrow(trainingSet3$problemSequences) #6557/95279
length(trainingSet3$problemGroups) #228
pg <- trainingSet3$problemGroups
pg_level <- sapply(strsplit(pg, ";"), function(x) sum(x != "")) - 1L
pg_df <- data.frame(ProblemGroup = pg, Level = pg_level, stringsAsFactors = FALSE)

norm <- function(x) gsub("\\s*;\\s*", ";", trimws(x))
pg_idx <- match(norm(trainingSet3$problemGroups), norm(trainingSet3$taxonomy))
sum(is.na(pg_idx))  # should be 0
pg_level_decipher <- trainingSet3$levels[pg_idx]
problem <- table(pg_level_decipher)
total   <- table(trainingSet3$levels)
out <- cbind(
  Problem = problem,
  Total   = total[names(problem)],
  Prop    = problem / total[names(problem)]
)
out

#  Problem Total  Proportion
#2       1     3 0.333333333
#3       4    15 0.266666667
#4      12   169 0.071005917
#5      30   402 0.074626866
#6      16   999 0.016016016
#7      13  1503 0.008649368
#8     152  4070 0.037346437

trainingSet3$problemGroups[pg_level_decipher == 2]
trainingSet3$problemGroups[pg_level_decipher == 3]
trainingSet3$problemGroups[pg_level_decipher == 5]
trainingSet3$problemGroups[pg_level_decipher == 6]

saveRDS(trainingSet3, file = "SILVA_LSU_138.2_trainingSet_ranked.rds")

#Without rank
trainingSet4 <- LearnTaxa(
  train    = ref_seqs,
  taxonomy = tax_strings,
  rank     = NULL,
  K        = NULL,       #NULL lets DECIPHER choose k‑mer size automatically
  maxIterations = 10,	 #Number of training iterations
  verbose  = TRUE
)
nrow(trainingSet4$problemSequences) #6545/95279
length(trainingSet4$problemGroups) #233
pg <- trainingSet4$problemGroups

norm <- function(x) gsub("\\s*;\\s*", ";", trimws(x))
pg_idx <- match(norm(trainingSet4$problemGroups), norm(trainingSet4$taxonomy))
sum(is.na(pg_idx))  # should be 0
pg_level_decipher <- trainingSet4$levels[pg_idx]
problem <- table(pg_level_decipher)
total   <- table(trainingSet4$levels)
out <- cbind(
  Problem = problem,
  Total   = total[names(problem)],
  Prop    = problem / total[names(problem)]
)
out

#  Problem Total        Prop
#2       1     3 0.333333333
#3       4    15 0.266666667
#4      14   169 0.082840237
#5      31   402 0.077114428
#6      16   999 0.016016016
#7      13  1503 0.008649368
#8     154  4070 0.037837838

trainingSet4$problemGroups[pg_level_decipher == 2]
trainingSet4$problemGroups[pg_level_decipher == 3]
trainingSet4$problemGroups[pg_level_decipher == 5]
trainingSet4$problemGroups[pg_level_decipher == 6]

saveRDS(trainingSet4, file = "SILVA_LSU_138.2_trainingSet_unranked.rds")
```
Run IDTAXA classification of the 18S and 28S ASVs
```R
library(DECIPHER)
library(Biostrings)

load("SILVA_SSU_r138_2_2024.RData")
SSU_trainingSet <- trainingSet
LSU_trainingSet <- readRDS("SILVA_LSU_138.2_trainingSet_ranked.rds")

SSU_ASVs <- readDNAStringSet("download_18122025/ASV_18S/ASVs.fasta")
LSU_ASVs <- readDNAStringSet("download_18122025/ASV_28S/ASVs.fasta")

SSU_tax_idtaxa <- IdTaxa(
  SSU_ASVs,
  SSU_trainingSet,
  strand = "both",

  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  #conservative - raise to 70–80 for fewer false positives
  verbose = TRUE 
)

# Convert to table
target_ranks <- c("root","domain","phylum","class","order","family","genus","species")

extract_tax <- function(x) {
  out <- setNames(rep(NA_character_, length(target_ranks)), target_ranks)
  if (!is.null(x$rank) && length(x$rank)) {
    idx <- match(tolower(x$rank), target_ranks)
    keep <- !is.na(idx)
    out[idx[keep]] <- x$taxon[keep]
  }
  out
}

SSU_tax_tab <- t(vapply(SSU_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
SSU_tax_tab <- as.data.frame(SSU_tax_tab, stringsAsFactors = FALSE)
rownames(SSU_tax_tab) <- names(SSU_ASVs)
write.table(SSU_tax_tab, file = "download_18122025/ASV_18S/SSU_tax_tab_corrected.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")

LSU_tax_idtaxa <- IdTaxa(
  LSU_ASVs,
  LSU_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  # conservative; raise to 70–80 if you want fewer false positives
  verbose = TRUE 
)

LSU_tax_tab <- t(vapply(LSU_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
LSU_tax_tab <- as.data.frame(LSU_tax_tab, stringsAsFactors = FALSE)
rownames(LSU_tax_tab) <- names(LSU_ASVs)

#Reverse ':' seperated taxonomic path conversion
extract_last <- function(x) {
  if (is.na(x)) {
    return(NA)
  }
  parts <- strsplit(x, ":")[[1]]
  return(parts[length(parts)])
}

# Apply the function to each element of each column
LSU_tax_tab[, 3:ncol(LSU_tax_tab)] <- lapply(LSU_tax_tab[, 3:ncol(LSU_tax_tab)], function(col) sapply(col, extract_last))

write.table(LSU_tax_tab, file = "download_18122025/ASV_28S/LSU_tax_tab_corrected.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")
```

### IDTAXA Custom <a name="19"></a>

Trial IDTAXA classification with the truncated dataset matching our primers (generated by QIIME2 (3.2)).

Export truncated SSU and LSU datasets to .fasta file and pull down for use in R with IDTAXA.
```bash
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path ~/db/SILVA/ssu_138.2/NEM-F-reference-amplicon-110.qza \
  --output-path ~/db/SILVA/ssu_138.2/silva-138.2-ssu-nr99-seqs/NEM-F

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path ~/db/SILVA/lsu_138.2/D2A-F-reference-amplicon-160.qza \
  --output-path ~/db/SILVA/lsu_138.2/silva-138.2-lsu-nr99-seqs/D2A-F
```
Train the IDTAXA classifier for the NEM-F SSU dataset:
```R
tax_table <- read.delim("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/down_20260116/silva-138.2-ssu-nr99-taxonomy.tsv",
                        sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
ref_seqs <- readDNAStringSet("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/down_20260116/NEM-F-dna-sequences.fasta")
taxid <- read.delim("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/down_20260116/tax_slv_ssu_138.2.txt", header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE,
                    col.names = c("hierarchy", "Index", "Rank", "tmp", "version"))
seq_id <- sub(" .*", "", names(ref_seqs))
m <- match(seq_id, tax_table$Feature.ID)
tax_table2 <- tax_table[m, ]
stopifnot(all(tax_table2$Feature.ID == seq_id))
all(seq_id == tax_table2$Feature.ID) #should return 'TRUE'

#Clean and prepare taxonomy strings
tax_clean <- str_replace_all(tax_table2$Taxon, "[a-z]__","")
tax_list <- str_split(tax_clean, ";\\s*")
tax_mat <- do.call(rbind, tax_list)
colnames(tax_mat) <- c("domain","kingdom","phylum","class", "order","family","genus","species")
tax_mat[tax_mat == ""] <- NA

tax_mat <- t(apply(tax_mat, 1, propagate_missing))

#all levels
tax_mat2 <- expand_tax_mat_to_paths(tax_mat)
tax_strings2 <- apply(tax_mat2, 1, function(r) {paste(c("Root", r[!is.na(r)]), collapse="; ")})
tax_strings2 <- trimws(tax_strings2)
tax_strings2 <- gsub("\\s*;\\s*", "; ", tax_strings2)

#Build rank taxid table from the taxonomy strings
level_to_rank <- c("root", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
spl <- strsplit(tax_strings2, ";\\s*")
edge_list <- do.call(rbind, lapply(spl, function(x) {
  x <- x[x != ""]
  paths <- vapply(seq_along(x), function(i) paste(x[1:i], collapse="; "), character(1))
  parents <- c(NA_character_, paths[-length(paths)])
  data.frame(
    path = paths,
    parent_path = parents,
    level = seq_along(paths) - 1L,           # Root = 0
    name = x,
    stringsAsFactors = FALSE
  )
}))
nodes <- unique(edge_list[, c("path", "parent_path", "level", "name")])
nodes <- nodes[order(nodes$level, nodes$path), ]
nodes$Index <- seq_len(nrow(nodes)) - 1L
path_to_index <- setNames(nodes$Index, nodes$path)
nodes$Parent <- ifelse(is.na(nodes$parent_path), -1L, path_to_index[nodes$parent_path])
nodes$Rank <- level_to_rank[pmin(nodes$level + 1L, length(level_to_rank))]
rank_table3 <- nodes[, c("Index", "name", "Parent", "level", "Rank")]
colnames(rank_table3) <- c("Index", "Name", "Parent", "Level", "Rank")

#Check for any remaining ambiguous names
ambiguous <- rank_table3[
  duplicated(rank_table3[, c("Name", "Level")]) |
  duplicated(rank_table3[, c("Name", "Level")], fromLast = TRUE),
  c("Name", "Level")
]
unique(ambiguous)

#Train the IDTAXA classifier
#With rank
SSU_trainingSet <- LearnTaxa(
  train    = ref_seqs,
  taxonomy = tax_strings2,
  rank     = rank_table3,
  K        = NULL,       #NULL lets DECIPHER choose k‑mer size automatically
  maxIterations = 10,	 #Number of training iterations
  verbose  = TRUE
)
nrow(SSU_trainingSet$problemSequences) #2/49441
length(SSU_trainingSet$problemGroups) #1
pg <- SSU_trainingSet$problemGroups
pg_level <- sapply(strsplit(pg, ";"), function(x) sum(x != "")) - 1L
pg_df <- data.frame(ProblemGroup = pg, Level = pg_level, stringsAsFactors = FALSE)

norm <- function(x) gsub("\\s*;\\s*", ";", trimws(x))
pg_idx <- match(norm(SSU_trainingSet$problemGroups), norm(SSU_trainingSet$taxonomy))
sum(is.na(pg_idx))  # should be 0
pg_level_decipher <- SSU_trainingSet$levels[pg_idx]
problem <- table(pg_level_decipher)
total   <- table(SSU_trainingSet$levels)
out <- cbind(
  Problem = problem,
  Total   = total[names(problem)],
  Prop    = problem / total[names(problem)]
)
out

#  Problem Total         Prop
#7       1  1433 0.0006978367

SSU_trainingSet$problemGroups[pg_level_decipher == 2]
SSU_trainingSet$problemGroups[pg_level_decipher == 3]
SSU_trainingSet$problemGroups[pg_level_decipher == 5]
SSU_trainingSet$problemGroups[pg_level_decipher == 6]

saveRDS(SSU_trainingSet, file = "SILVA_SSU_138.2_NEM-F_trainingSet_ranked.rds")
```
Train the IDTAXA classifier for the D2A-F LSU dataset:
```R
#Import data
tax_table <- read.delim("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/down_20260116/silva-138.2-lsu-nr99-taxonomy.tsv",
                        sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
ref_seqs <- readDNAStringSet("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/down_20260116/D2A-F-dna-sequences.fasta")
taxid <- read.delim("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/down_20260116/tax_slv_lsu_138.2.txt", header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE,
                    col.names = c("hierarchy", "Index", "Rank", "tmp", "version"))
seq_id <- sub(" .*", "", names(ref_seqs))
m <- match(seq_id, tax_table$Feature.ID)
tax_table2 <- tax_table[m, ]
stopifnot(all(tax_table2$Feature.ID == seq_id))
all(seq_id == tax_table2$Feature.ID) #should return 'TRUE'

#Clean and prepare taxonomy strings
tax_clean <- str_replace_all(tax_table2$Taxon, "[a-z]__","")
tax_list <- str_split(tax_clean, ";\\s*")
tax_mat <- do.call(rbind, tax_list)
colnames(tax_mat) <- c("domain","kingdom","phylum","class", "order","family","genus","species")
tax_mat[tax_mat == ""] <- NA

tax_mat <- t(apply(tax_mat, 1, propagate_missing))

#all levels
tax_mat2 <- expand_tax_mat_to_paths(tax_mat)
tax_strings2 <- apply(tax_mat2, 1, function(r) {paste(c("Root", r[!is.na(r)]), collapse="; ")})
tax_strings2 <- trimws(tax_strings2)
tax_strings2 <- gsub("\\s*;\\s*", "; ", tax_strings2)

#Build rank taxid table from the taxonomy strings
level_to_rank <- c("root", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
spl <- strsplit(tax_strings2, ";\\s*")
edge_list <- do.call(rbind, lapply(spl, function(x) {
  x <- x[x != ""]
  paths <- vapply(seq_along(x), function(i) paste(x[1:i], collapse="; "), character(1))
  parents <- c(NA_character_, paths[-length(paths)])
  data.frame(
    path = paths,
    parent_path = parents,
    level = seq_along(paths) - 1L,           # Root = 0
    name = x,
    stringsAsFactors = FALSE
  )
}))
nodes <- unique(edge_list[, c("path", "parent_path", "level", "name")])
nodes <- nodes[order(nodes$level, nodes$path), ]
nodes$Index <- seq_len(nrow(nodes)) - 1L
path_to_index <- setNames(nodes$Index, nodes$path)
nodes$Parent <- ifelse(is.na(nodes$parent_path), -1L, path_to_index[nodes$parent_path])
nodes$Rank <- level_to_rank[pmin(nodes$level + 1L, length(level_to_rank))]
rank_table3 <- nodes[, c("Index", "name", "Parent", "level", "Rank")]
colnames(rank_table3) <- c("Index", "Name", "Parent", "Level", "Rank")

#Check for any remaining ambiguous names
ambiguous <- rank_table3[
  duplicated(rank_table3[, c("Name", "Level")]) |
  duplicated(rank_table3[, c("Name", "Level")], fromLast = TRUE),
  c("Name", "Level")
]
unique(ambiguous)

#Train the IDTAXA classifier
#With rank
LSU_trainingSet <- LearnTaxa(
  train    = ref_seqs,
  taxonomy = tax_strings2,
  rank     = rank_table3,
  K        = NULL,       #NULL lets DECIPHER choose k‑mer size automatically
  maxIterations = 10,	 #Number of training iterations
  verbose  = TRUE
)
nrow(LSU_trainingSet$problemSequences) #7/4662
length(LSU_trainingSet$problemGroups) #5
pg <- LSU_trainingSet$problemGroups
pg_level <- sapply(strsplit(pg, ";"), function(x) sum(x != "")) - 1L
pg_df <- data.frame(ProblemGroup = pg, Level = pg_level, stringsAsFactors = FALSE)

norm <- function(x) gsub("\\s*;\\s*", ";", trimws(x))
pg_idx <- match(norm(LSU_trainingSet$problemGroups), norm(LSU_trainingSet$taxonomy))
sum(is.na(pg_idx))  # should be 0
pg_level_decipher <- LSU_trainingSet$levels[pg_idx]
problem <- table(pg_level_decipher)
total   <- table(LSU_trainingSet$levels)
out <- cbind(
  Problem = problem,
  Total   = total[names(problem)],
  Prop    = problem / total[names(problem)]
)
out

#        Problem  Total       Prop
#5       5        78          0.06410256

LSU_trainingSet$problemGroups[pg_level_decipher == 2]
LSU_trainingSet$problemGroups[pg_level_decipher == 3]
LSU_trainingSet$problemGroups[pg_level_decipher == 5]
LSU_trainingSet$problemGroups[pg_level_decipher == 6]

saveRDS(LSU_trainingSet, file = "SILVA_LSU_138.2_D2A-F_trainingSet_ranked.rds")
```
Run IDTAXA classification of the 18S and 28S ASVs
```R
#Import
SSU_trainingSet <- readRDS("SILVA_SSU_138.2_NEM-F_trainingSet_ranked.rds")
LSU_trainingSet <- readRDS("SILVA_LSU_138.2_D2A-F_trainingSet_ranked.rds")

SSU_ASVs <- readDNAStringSet("download_18122025/ASV_18S/ASVs.fasta")
LSU_ASVs <- readDNAStringSet("download_18122025/ASV_28S/ASVs.fasta")

#SSU
SSU_tax_idtaxa <- IdTaxa(
  SSU_ASVs,
  SSU_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  #conservative - raise to 70–80 for fewer false positives
  verbose = TRUE 
)

target_ranks <- c("root","domain","phylum","class","order","family","genus","species")
SSU_tax_tab <- t(vapply(SSU_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
SSU_tax_tab <- as.data.frame(SSU_tax_tab, stringsAsFactors = FALSE)
rownames(SSU_tax_tab) <- names(SSU_ASVs)

#LSU
LSU_tax_idtaxa <- IdTaxa(
  LSU_ASVs,
  LSU_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  # conservative; raise to 70–80 if you want fewer false positives
  verbose = TRUE 
)

LSU_tax_tab <- t(vapply(LSU_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
LSU_tax_tab <- as.data.frame(LSU_tax_tab, stringsAsFactors = FALSE)
rownames(LSU_tax_tab) <- names(LSU_ASVs)

#Export
#Reverse ':' seperated taxonomic path conversion
SSU_tax_tab[, 3:ncol(SSU_tax_tab)] <- lapply(SSU_tax_tab[, 3:ncol(SSU_tax_tab)], function(col) sapply(col, extract_last))
LSU_tax_tab[, 3:ncol(LSU_tax_tab)] <- lapply(LSU_tax_tab[, 3:ncol(LSU_tax_tab)], function(col) sapply(col, extract_last))

write.table(SSU_tax_tab, file = "download_18122025/ASV_18S/NEM-F-SSU_tax_tab_corrected.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")
write.table(LSU_tax_tab, file = "download_18122025/ASV_28S/D2A-F-LSU_tax_tab_corrected.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")
```
## Taxonomy Assignment - BLASTN <a name="5"></a>
### BLASTN <a name="14"></a>

BLAST is alignment based, limited by 'best hit' interpretation, and prone to misidentifying short or conserved seqeunces. For eDNA/metabarcoding, BLAST is often too literal — it finds the closest sequence, even if it’s wrong. However, the NCBI nt database is far larger than any dedicated database, including SILVA. BLAST against the nt database is the approach used previously by Veronika.

#### NCBI nt <a name="20"></a>
BLASTN vs NCBI nt database:

```bash
for ASV in $(find /data/users/theaven/nematode_project/asvs -name 'ASVs.fasta' -type f); do
	Task=blast
	Database=/data/blobtoolkit/nt/nt
	Max_target=10
	OutPrefix=$(dirname $ASV | rev | cut -d '/' -f1 | rev)
	OutDir="$(dirname $ASV)"/"$Task"
	mkdir -p $OutDir
	ExpectedOutput="$OutDir"/${OutPrefix}.vs."$(basename $Database)".mts"$Max_target".hsp1.1e25.megablast.out

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 5s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_blastn.sh "$ASV" "$Database" "$OutDir" "$OutPrefix" "$Max_target")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done

#Top hits extracted to:
/data/users/theaven/nematode_project/asvs/ASV_18S/blast/ASV_18S.vs.nt.mts1.hsp1.1e25.megablast.out
/data/users/theaven/nematode_project/asvs/ASV_28S/blast/ASV_28S.vs.nt.mts1.hsp1.1e25.megablast.out

module load anaconda3
conda activate taxonkit
wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz
mkdir -p $HOME/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit

cut -f 2 /data/users/theaven/nematode_project/asvs/ASV_18S/blast/ASV_18S.vs.nt.mts1.hsp1.1e25.megablast.out > taxonomy_ids.txt
taxonkit lineage taxonomy_ids.txt | taxonkit reformat -F -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" > taxonomic_paths.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR{gsub(/\r/,"",$0);tax[$1]=$2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8;next}{gsub(/\r/,"",$0)}FNR==1{print $0,"k","p","c","o","f","g","s";next}{tid=$2;if(tid~/;/)print $0,"","","","","","","";else if(tid in tax)print $0,tax[tid];else print $0,"","","","","","",""}' taxonomic_paths.tsv /data/users/theaven/nematode_project/asvs/ASV_18S/blast/ASV_18S.vs.nt.mts1.hsp1.1e25.megablast.out > /data/users/theaven/nematode_project/asvs/ASV_18S/blast/ASV_18S.megablast.with_taxonomy.tsv

cut -f 2 /data/users/theaven/nematode_project/asvs/ASV_28S/blast/ASV_28S.vs.nt.mts1.hsp1.1e25.megablast.out > taxonomy_ids.txt
taxonkit lineage taxonomy_ids.txt | taxonkit reformat -F -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" > taxonomic_paths.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR{gsub(/\r/,"",$0);tax[$1]=$2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8;next}{gsub(/\r/,"",$0)}FNR==1{print $0,"k","p","c","o","f","g","s";next}{tid=$2;if(tid~/;/)print $0,"","","","","","","";else if(tid in tax)print $0,tax[tid];else print $0,"","","","","","",""}' taxonomic_paths.tsv /data/users/theaven/nematode_project/asvs/ASV_28S/blast/ASV_28S.vs.nt.mts1.hsp1.1e25.megablast.out > /data/users/theaven/nematode_project/asvs/ASV_28S/blast/ASV_28S.megablast.with_taxonomy.tsv
```
#### SILVA SSU and LSU - BLAST locally / ACT lookup <a name="21"></a>

SILVA ACT - SILVA’s ACT service allows you to upload your own sequences and have them aligned, taxonomically classified, and placed into a reference phylogenetic tree. This is the official and recommended way to classify your 18S/28S sequences using SILVA. SINA aligner (SILVA’s standalone alignment tool) is an alignment based classifier, similar to mothur / ARB.

ASV seqeunces were submitted to the ACT tool for annotation. For comparison ASVs were also BLASTN searched against the SILVA database.
```bash
#Convert SILVA to DNA
module load anaconda3
conda activate seqkit
seqkit seq --remove-gaps --rna2dna ~/db/SILVA/ssu_138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta > ~/db/SILVA/ssu_138.2/SILVA_138.2_SSURef_NR99_tax_silva_dna.fasta
seqkit seq --remove-gaps --rna2dna ~/db/SILVA/lsu_138.2/SILVA_138.2_LSURef_NR99_tax_silva.fasta > ~/db/SILVA/lsu_138.2/SILVA_138.2_LSURef_NR99_tax_silva_dna.fasta

#Build BLAST databases for SILVA SSU and LSU
conda activate blast
cd ~/db/blastn
makeblastdb -in ~/db/SILVA/ssu_138.2/SILVA_138.2_SSURef_NR99_tax_silva_dna.fasta -dbtype nucl -parse_seqids -out silva_ssu
makeblastdb -in ~/db/SILVA/lsu_138.2/SILVA_138.2_LSURef_NR99_tax_silva_dna.fasta -dbtype nucl -parse_seqids -out silva_lsu

#Run BLASTN
for ASV in /data/users/theaven/nematode_project/asvs/ASV_18S/ASVs.fasta; do
	Task=blast
	Database=~/db/blastn/silva_ssu
	Max_target=1
	OutPrefix=$(dirname $ASV | rev | cut -d '/' -f1 | rev)
	OutDir="$(dirname $ASV)"/"$Task"
	mkdir -p $OutDir
	ExpectedOutput="$OutDir"/${OutPrefix}.vs."$(basename $Database)".mts"$Max_target".hsp1.1e25.megablast.out

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 5s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_blastn.sh "$ASV" "$Database" "$OutDir" "$OutPrefix" "$Max_target")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done

for ASV in /data/users/theaven/nematode_project/asvs/ASV_28S/ASVs.fasta; do
	Task=blast
	Database=~/db/blastn/silva_lsu
	Max_target=1
	OutPrefix=$(dirname $ASV | rev | cut -d '/' -f1 | rev)
	OutDir="$(dirname $ASV)"/"$Task"
	mkdir -p $OutDir
	ExpectedOutput="$OutDir"/${OutPrefix}.vs."$(basename $Database)".mts"$Max_target".hsp1.1e25.megablast.out

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 9 ]; do
		sleep 5s
		printf "."
		Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	done

	if [ ! -s "$ExpectedOutput" ]; then
		jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_blastn.sh "$ASV" "$Database" "$OutDir" "$OutPrefix" "$Max_target")
		printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
	else
		echo "For $ID found: $ExpectedOutput" 
	fi
done
```

## Taxonomy Assignment - BOLD <a name="22"></a>

"BOLD enables researchers, students and citizen scientists to submit DNA barcode sequences for identification against a reference library that includes sequences from more than 1.2 million species." The BOLD search method used k-mer/distance-based matching, its own alignment and similarity scoring (optimized for COI), and BIN system (Barcode Index Number clusters—computed with algorithms not based on BLAST). The BOLD database is extensive for COI but far less deep for 18S seqeunces, 28S depth is lower still and 28S is not incuded in any of the searchable datasets.

The 18S ASV sequences were searched against the "FUNGI LIBRARY (PUBLIC)" and "ANIMAL SECONDARY MARKERS (PUBLIC)" datasets with "Exhaustive Search" operating mode [16th-Jan-2026] using the BOLDSYSTEMS web portal https://id.boldsystems.org/

## Nematode Curated Databases <a name="23"></a>

Several curated nematode specific barcode datasets have been published:

NemaTaxa - Built explicitly from sequences that match the NF1 (5′-GCCTCCCTCGCGCCATCAGGGTGGTGCATGGCCGTTCTTAGTT-3′) (forward) / 18Sr2b 18Sr2b (5′-GCCTTGCCAGCCCGCTCAGTACAAAGGGCAGGGACGTAAT-3’) (reverse) primer pair (the standard soil nematode metabarcoding marker). Source: Nematode 18S sequences pulled from NCBI + SILVA v132. Curation: Manual trimming of taxonomy to classical Linnaean ranks (kingdom → genus). Missing ranks filled manually; inconsistent strings cleaned. Good for the the NF1 / 18Sr2b primer pair soil/agricultural metabarcoding. Very clean taxonomy strings suitable for QIIME2.

NemaBase - “nematode-optimized” subset of 18S rRNA reference sequences. Source: all nematode 18S sequences from SILVA (v111 and v138). Two versions: full length and NF1 /18Sr2b primer pair. Curation: cleaning up taxonomy (standardizing taxonomic ranks, updating according to accepted nematode classification, removing suspicious or erroneous entries) to produce a cleaner, well-annotated dataset. Best for soil, freshwater or terrestrial nematode communities where sequenced nematodes are common in SILVA and common metabarcoding primers are used.

Nemabiome - Charrier et al. 2024 rRNA-cistron 18S DB - Four databases; 18S, 28S, ITS-1, and ITS-1-5.8S-ITS-2. Curation: Built with the markerDB pipeline (systematic mining and filtering of NCBI GenBank) and deduplicated to full-length representatives. Has maximum species coverage across Nematoda. Best for wildlife parasites, marine/free-living nematodes, or very broad community surveys where phylum-wide coverage matters.


### IDTAXA - nemabiome <a name="24"></a>

It seem that the nemabiome dataset is most relavent as we have not used the NF1 / 18Sr2b primer pair, we have 28S seqeunces, and we are looking at wildlife parasites. Also as classification has already been performed with the SILVA database a nematode dataset pulled from NCBI GenBank is potentially more addative. IDTAXA files were downloaded from https://www.nemabiome.ca/ [Accessed: 16th-Jan-2026].
18S = 2,645 non-redundant sequences, 16 orders, 178 families, 514 genera, 1,391 species
28S = 254 non-redundant sequences, 4 orders, 44 families, 117 genera, 204 species

Train classifiers from the nemabiome datasets:
```R
library(DECIPHER)
library(Biostrings)
library(stringr)

setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")


#SSU
ref_seqs <- readDNAStringSet("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/nemabiome_idtaxa_18S.fasta")
ref_tax <- names(ref_seqs)
taxid <- read.delim("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/nemabiome_idtaxa_18S.tax", sep = "\t", header = TRUE)
SSU_nema_trainingSet <- LearnTaxa(
  train    = ref_seqs,
  taxonomy = ref_tax,
  rank     = taxid,
  K        = NULL,       #NULL lets DECIPHER choose k‑mer size automatically
  maxIterations = 10,	 #Number of training iterations
  verbose  = TRUE
)
saveRDS(SSU_nema_trainingSet, file = "nemabiome_18S_trainingSet_unranked.rds")

#LSU
ref_seqs <- readDNAStringSet("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/nemabiome_idtaxa_28s.fasta")
ref_tax <- names(ref_seqs)
taxid <- read.delim("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/nemabiome_idtaxa_28s.tax", sep = "\t", header = TRUE)
LSU_nema_trainingSet <- LearnTaxa(
  train    = ref_seqs,
  taxonomy = ref_tax,
  rank     = taxid,
  K        = NULL,       #NULL lets DECIPHER choose k‑mer size automatically
  maxIterations = 10,	 #Number of training iterations
  verbose  = TRUE
)
saveRDS(LSU_nema_trainingSet, file = "nemabiome_28S_trainingSet_unranked.rds")
```
Taxonomic classification with nemabiome datasets via IDTAXA:
```R
SSU_ASVs <- readDNAStringSet("download_18122025/ASV_18S/ASVs.fasta")
LSU_ASVs <- readDNAStringSet("download_18122025/ASV_28S/ASVs.fasta")

#SSU
SSU_tax_idtaxa <- IdTaxa(
  SSU_ASVs,
  SSU_nema_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  #conservative - raise to 70–80 for fewer false positives
  verbose = TRUE 
)

SSU_tax_tab <- t(vapply(SSU_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
SSU_tax_tab <- as.data.frame(SSU_tax_tab, stringsAsFactors = FALSE)
rownames(SSU_tax_tab) <- names(SSU_ASVs)
write.table(SSU_tax_tab, file = "download_18122025/ASV_18S/nemabiome_SSU_tax_tab_corrected.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")


#LSU
LSU_tax_idtaxa <- IdTaxa(
  LSU_ASVs,
  LSU_nema_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  #conservative - raise to 70–80 for fewer false positives
  verbose = TRUE 
)

LSU_tax_tab <- t(vapply(LSU_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
LSU_tax_tab <- as.data.frame(LSU_tax_tab, stringsAsFactors = FALSE)
rownames(LSU_tax_tab) <- names(LSU_ASVs)
write.table(LSU_tax_tab, file = "download_18122025/ASV_28S/nemabiome_LSU_tax_tab_corrected.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")
```

### IDTAXA - nemabase <a name="24"></a>

The nemabase-18S dataset was downloaded from https://github.com/WormsEtAl/18SNemaBase [Accessed 16th-Jan-2026]

Train classifier from the nemabase datasets:
```R
library(DECIPHER)
library(Biostrings)
library(stringr)

setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")

ref_seqs <- readDNAStringSet("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R/18S_NemaBase.fasta")
tax_strings <- sub("^[^_]+_", "", names(ref_seqs))
tax_strings <- ifelse(grepl("^Root;", tax_strings), tax_strings, paste0("Root;", tax_strings))
tax_strings <- gsub("\\s*;\\s*", ";", tax_strings)

expand_one <- function(tx) {
  parts <- strsplit(tx, ";", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]

  expanded <- vapply(
    seq_along(parts),
    function(i) paste(parts[1:i], collapse = ":"),
    character(1)
  )

  paste(expanded, collapse = ";")
}

tax_strings_expanded <- vapply(tax_strings, expand_one, character(1))


#Build rank taxid table from the taxonomy strings
level_to_rank <- c("root", "domain", "kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "superorder", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species")
spl <- lapply(strsplit(tax_strings_expanded, ";"), trimws)
edge_list <- do.call(rbind, lapply(spl, function(x) {
  x <- x[nzchar(x)]
  paths <- vapply(seq_along(x), function(i) paste(x[1:i], collapse=";"), character(1))
  parents <- c(NA_character_, paths[-length(paths)])
  data.frame(
    path = paths,
    parent_path = parents,
    level = seq_along(paths) - 1L,  # Root = 0
    name = x,
    stringsAsFactors = FALSE
  )
}))

nodes <- unique(edge_list[, c("path", "parent_path", "level", "name")])
nodes <- nodes[order(nodes$level, nodes$path), ]
nodes$Index <- seq_len(nrow(nodes)) - 1L
path_to_index <- setNames(nodes$Index, nodes$path)
nodes$Parent <- ifelse(is.na(nodes$parent_path), -1L, unname(path_to_index[nodes$parent_path]))
nodes$Rank <- level_to_rank[pmin(nodes$level + 1L, length(level_to_rank))]
rank_table <- nodes[, c("Index", "name", "Parent", "level", "Rank")]
colnames(rank_table) <- c("Index", "Name", "Parent", "Level", "Rank")

#Check for any remaining ambiguous names
ambiguous <- rank_table[
  duplicated(rank_table[, c("Name", "Level")]) |
  duplicated(rank_table[, c("Name", "Level")], fromLast = TRUE),
  c("Name", "Level")
]
unique(ambiguous)

SSU_base_trainingSet <- LearnTaxa(
  train    = ref_seqs,
  taxonomy = tax_strings_expanded,
  rank     = rank_table,
  K        = NULL,       #NULL lets DECIPHER choose k‑mer size automatically
  maxIterations = 10,	 #Number of training iterations
  verbose  = TRUE
)
saveRDS(SSU_base_trainingSet, file = "nemabase_18S_trainingSet_ranked.rds")
```
Taxonomic classification with nemabase datasets via IDTAXA:
```R
SSU_ASVs <- readDNAStringSet("download_18122025/ASV_18S/ASVs.fasta")

SSU_tax_idtaxa <- IdTaxa(
  SSU_ASVs,
  SSU_base_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  #conservative - raise to 70–80 for fewer false positives
  verbose = TRUE 
)

target_ranks <- c("root", "domain", "kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "superorder", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species")
SSU_tax_tab <- t(vapply(SSU_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
SSU_tax_tab <- as.data.frame(SSU_tax_tab, stringsAsFactors = FALSE)
rownames(SSU_tax_tab) <- names(SSU_ASVs)
SSU_tax_tab[, 3:ncol(SSU_tax_tab)] <- lapply(SSU_tax_tab[, 3:ncol(SSU_tax_tab)], function(col) sapply(col, extract_last))
write.table(SSU_tax_tab, file = "download_18122025/ASV_18S/nemabase_SSU_tax_tab_corrected.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")
```
## Parse classification outputs
```bash
sed -e 's/";"/\t/g' -e 's/^"//' -e 's/"$//' -e 's/";/\t/' /data/users/theaven/nematode_project/asvs/ASV_18S/ACT/ACT-SILVA-18S.csv | tr -d '"' | sed -e 's/;/\t/' > /data/users/theaven/nematode_project/asvs/ASV_18S/ACT/ACT-SILVA-18S.tsv
sed -e 's/";"/\t/g' -e 's/^"//' -e 's/"$//' -e 's/";/\t/' /data/users/theaven/nematode_project/asvs/ASV_28S/ACT/ACT-SILVA-28S.csv | tr -d '"' | sed -e 's/;/\t/' > /data/users/theaven/nematode_project/asvs/ASV_28S/ACT/ACT-SILVA-28S.tsv

module load apptainer/1.4.1-gcc-13.3.0-3coysxn
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/unibz/merge_asv_classifications.py \
  --input "QIIME - full length:/data/users/theaven/nematode_project/asvs/ASV_18S/qiime_fullength_taxonomy/taxonomy.tsv:1" \
  --input "QIIME - custom:/data/users/theaven/nematode_project/asvs/ASV_18S/qiime_NEM-F_taxonomy/taxonomy.tsv:1" \
  --input "IDTAXA - full length:/data/users/theaven/nematode_project/asvs/ASV_18S/idtaxa/SSU_tax_tab_corrected.tsv:1" \
  --input "IDTAXA - custom:/data/users/theaven/nematode_project/asvs/ASV_18S/idtaxa/NEM-F-SSU_tax_tab_corrected.tsv:1" \
  --input "BLAST -nt:/data/users/theaven/nematode_project/asvs/ASV_18S/blast/ASV_18S.megablast.with_taxonomy.tsv:1" \
  --input "BLAST - SILVA:/data/users/theaven/nematode_project/asvs/ASV_18S/blast/ASV_18S.vs.silva_ssu.mts1.hsp1.1e25.megablast.out:1" \
  --input "ACT - SILVA:/data/users/theaven/nematode_project/asvs/ASV_18S/ACT/ACT-SILVA-18S.tsv:3" \
  --input "BOLD:/data/users/theaven/nematode_project/asvs/ASV_18S/BOLD/BOLD-18S.tsv:1" \
  --input "NEMABIOME:/data/users/theaven/nematode_project/asvs/ASV_18S/idtaxa/nemabiome_SSU_tax_tab_corrected.tsv:1" \
  --input "NEMABASE:/data/users/theaven/nematode_project/asvs/ASV_18S/idtaxa/nemabase_SSU_tax_tab_corrected.tsv:1" \
  -o merged_ASV_classifications_18S.tsv

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/unibz/merge_asv_classifications.py \
  --input "QIIME - full length:/data/users/theaven/nematode_project/asvs/ASV_28S/qiime_fullength_taxonomy/taxonomy.tsv:1" \
  --input "QIIME - custom:/data/users/theaven/nematode_project/asvs/ASV_28S/qiime_D2A-F_taxonomy/taxonomy.tsv:1" \
  --input "IDTAXA - full length:/data/users/theaven/nematode_project/asvs/ASV_28S/idtaxa/LSU_tax_tab_corrected.tsv:1" \
  --input "IDTAXA - custom:/data/users/theaven/nematode_project/asvs/ASV_28S/idtaxa/D2A-F-LSU_tax_tab_corrected.tsv:1" \
  --input "BLAST -nt:/data/users/theaven/nematode_project/asvs/ASV_28S/blast/ASV_28S.vs.nt.mts1.hsp1.1e25.megablast.out:1" \
  --input "BLAST - SILVA:/data/users/theaven/nematode_project/asvs/ASV_28S/blast/ASV_28S.vs.silva_lsu.mts1.hsp1.1e25.megablast.out:1" \
  --input "ACT - SILVA:/data/users/theaven/nematode_project/asvs/ASV_28S/ACT/ACT-SILVA-28S.tsv:3" \
  --input "NEMABIOME:/data/users/theaven/nematode_project/asvs/ASV_28S/idtaxa/nemabiome_LSU_tax_tab_corrected.tsv:1" \
  -o merged_ASV_classifications_28S.tsv
```

## Relate my ASVs to Veronika/Maja's reads <a name="26"></a>
```bash
#Combine read data, retaining sample info
for file in $(ls /data/users/theaven/nematode_project/Maja_data/*/*.fasta); do
Sample=$(basename $file | cut -d '_' -f1,2,3)
sed "/^>/ s/$/|$Sample/" $file >> merged_Maja_reads.fasta
done

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/substring_match.py \
  --short asvs/ASV_18S/ASVs.fasta \
  --long merged_Maja_reads.fasta \
  -o ASV18S_in_merged.tsv

grep -c '^>' merged_Maja_reads.fasta #3022183

#Get unique reads
module load apptainer/1.4.1-gcc-13.3.0-3coysxn
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/dedup_fasta.py merged_Maja_reads.fasta
#Input reads: 3022183
#Unique reads: 981511
#Wrote: merged_Maja_reads.fasta.unique.fasta
#Wrote: merged_Maja_reads.fasta.unique.map.tsv

grep -c '^>' asvs/ASV_18S/ASVs.fasta #629

#Find ASVs in the reads
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/substring_match.py \
  --short asvs/ASV_18S/ASVs.fasta \
  --long merged_Maja_reads.fasta.unique.fasta \
  -o ASV18S_in_merged_uniques.tsv

#How many ASV matches do reads have?:
awk -F'\t' 'NR>1{seen[$3 SUBSEP $1]=1} END{for(k in seen){split(k,a,SUBSEP); n[a[1]]++} for(l in n) d[n[l]]++ ; for(c in d) print c"\t"d[c]}' ASV18S_in_merged_uniques.tsv | sort -n
#1       386875
#2       4414
#3       353
#4       84
#5       30

awk -F'\t' -v OFS='\t' 'NR==1{next}{long=$3;asv=$1;key=long SUBSEP asv;if(!(seen[key]++))list[long]=list[long](list[long]?";":"")asv}END{for(long in list){n=split(list[long],a,";");if(n>1)print long,n,list[long]}}' ASV18S_in_merged_uniques.tsv | sort -k2,2nr




























#Load in Maja's BOLD annotations
module load anaconda3
conda activate xlsx2csv
for file in $(ls Maja_data/*/*.xlsx); do
	out=$(echo "$file" | sed 's@.xlsx@.tsv@g')
xlsx2csv "$file" | sed 's/,/\t/g' > "$out"
done

gawk -v FS='[[:space:]]+' -v OFS='\t' '
function na8(   i,s){ s="NA"; for(i=2;i<=8;i++) s=s OFS "NA"; return s }

# ---------- index all Maja files (everything except the last input file) ----------
ARGIND < ARGC-1 {
  fn = FILENAME
  sub(/^.*\//, "", fn)                 # basename

  sample = fn
  if (sub(/_trimmed.*/, "", sample) == 0) sub(/_unmapped.*/, "", sample)
  # sample is now filename prefix, e.g. NEM_ST2_10_1 or NEM_RV086_1

  id = $1
  if (id == "" || NF < 8) next

  map[sample SUBSEP id] = $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $1
  next
}

# ---------- process ASV file (the last input file) ----------
ARGIND == ARGC-1 {
  if (FNR == 1) next

  # Prefer tab-splitting to preserve "3rd tab column" even if it contains spaces
  nt = split($0, f, "\t")
  if (nt < 3) nt = split($0, f, /[[:space:]]+/)  # fallback if no tabs

  short_id = f[1]

  # Find rep=... anywhere on the full line (works regardless of which field it sits in)
  if (!match($0, /rep=[^[:space:]]+/, m)) next
  rep = m[0]
  sub(/^rep=/, "", rep)

  split(rep, a, /\|/)
  read_id = a[1]
  sample0 = a[2]          # e.g. NEM_ST2_10  (missing _1) OR NEM_RV086_1 (complete)

  # Look up: try exact sample, then try sample+"_1" (for ST1/ST2 cases)
  key = sample0 SUBSEP read_id
  if (!(key in map)) {
    key2 = (sample0 "_1") SUBSEP read_id
    if (key2 in map) key = key2
  }

  if (key in map) {
    print short_id, "Maja_BOLD", "NA", map[key], sample0
  } else {
    print short_id, "Maja_BOLD", "NA", na8(), sample0
  }
}
' Maja_data/*/*_unmapped_sorted_identification_result*.tsv ASV18S_in_merged_uniques.tsv \
> ASV_vs_reads_BOLD.tsv

cut -f3 ASV_vs_reads_BOLD.tsv | sort -u

awk -F'\t' '$4 != "no-match"' ASV_vs_reads_BOLD.tsv | awk -F'\t' '$4 == "NA"' | cut -f12 | sort -u
awk -F'\t' '$4 != "no-match"' ASV_vs_reads_BOLD.tsv > ASV_vs_reads_BOLD_informative.tsv

```

## Analyse <a name="27"></a>

Importing into R
```R
install.packages(c(
  "readxl",
  "tidyverse",
  "janitor",
  "vegan",
  "irr",
  "writexl",
  "dplyr", 
  "tidyr", 
  "purrr",
  "tibble",
  "stringr",
  "ggrepel"
))

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("phyloseq")
}

library(readxl)
library(tidyverse)
library(janitor)
library(vegan)
library(irr)
library(writexl)
library(dplyr)
library(purrr)
library(tidyr)
library(phyloseq)
library(tibble)
library(stringr)
library(ggrepel)

setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")
set.seed(1)

asv_file <- "merged_ASV_classifications_6.xlsx"

#ASV pipeline inputs:
asv_28S_long     <- read_xlsx(asv_file, sheet = "28S",                    na = c("", "NA")) |> clean_names()
asv_28S_long <- asv_28S_long |> mutate(sample_id = sample |> str_remove("^28S_") |> str_remove("^D2A_")|> str_to_lower())
asv_28S_counts   <- read_xlsx(asv_file, sheet = "28S_ASVs_to_samples",    na = c("", "NA")) |> clean_names()
asv_28S_counts <- asv_28S_counts |> rename_with( ~ .x |> str_remove("^x28s_d2a_") |> str_remove("^x28s_"), -feature_id )
manual_28S       <- read_xlsx(asv_file, sheet = "28S_tall_manual_only",   na = c("", "NA")) |> clean_names()

asv_18S_long     <- read_xlsx(asv_file, sheet = "18S",                    na = c("", "NA")) |> clean_names()
asv_18S_long <- asv_18S_long |> mutate(sample_id = sample |> str_remove("^18S_") |> str_remove("^NEM_")|> str_to_lower())
asv_18S_counts   <- read_xlsx(asv_file, sheet = "18S_ASVs_to_samples",    na = c("", "NA")) |> clean_names()
asv_18S_counts <- asv_18S_counts |> rename_with( ~ .x |> str_remove("^x18s_nem_") |> str_remove("^NEM_"), -feature_id )
manual_18S       <- read_xlsx(asv_file, sheet = "18S_tall_manual_only",   na = c("", "NA")) |> clean_names()


name_fixes <- c(
  "\\bDiplogasteridae\\b" = "Diplogastridae"
)

for (nm in c("asv_28S_long", "asv_28S_counts", "manual_28S",
             "asv_18S_long", "asv_18S_counts", "manual_18S")) {
  assign(nm,
         get(nm) |>
           mutate(across(where(is.character), \(x) str_replace_all(x, name_fixes))))
}

#Read pipeline inputs
veg <- "vegan_input_all.xlsx" 
sample_metadata <- read_xlsx(veg, "sample_metadata") |>
  mutate(yr    = coalesce(suppressWarnings(as.integer(year)),
                          2000L + as.integer(str_extract(group, "\\d{2}$"))),
         month = as.integer(str_match(group, "^[A-Z]_(\\d{2})_\\d{2}$")[, 2]),
         region = if_else(location %in% c("Abetone","Croatia","Romania","Russia"),
                          "non-alpine", "alpine"))


sample_metadata <- read_xlsx(veg, "sample_metadata") |>
  clean_names() |>
  mutate(
    sample_id = str_to_lower(sample_id),
    yr = coalesce(
      suppressWarnings(as.integer(year)),
      2000L + as.integer(str_extract(group, "\\d{2}$"))
    ),
    month = as.integer(
      str_match(group, "^[A-Z]_(\\d{2})_\\d{2}$")[, 2]
    ),
    country = case_when(
      location %in% c("Innichen", "Tiers", "Abetone", "Prad", "Italy") ~ "Italy",
      location %in% c("Kals", "Austria") ~ "Austria",
      location %in% c("Lika", "Croatia") ~ "Croatia",
      location %in% c("Cluj", "Romania") ~ "Romania",
      location %in% c("Moscow", "Russia") ~ "Russia",
      TRUE ~ NA_character_
    ),
    locality = case_when(
      location %in% c("Innichen", "Tiers", "Abetone", "Prad") ~ location,
      location %in% c("Lika", "Croatia") ~ "Lika",
      location %in% c("Cluj", "Romania") ~ "Cluj",
      location %in% c("Moscow", "Russia") ~ "Moscow",
      location == "Kals" ~ "Kals",
      TRUE ~ NA_character_
    ),
    region = if_else(
      country %in% c("Croatia", "Romania", "Russia"),
      "non-alpine",
      "alpine"
    )
  )



nri <- "Nematode_reads_info.xlsx"
clean_tax <- function(x) x |> str_squish() |> str_replace_all("[ _]+", "_") |>
  str_replace("sp\\._?O717", "sp._O717")
reads <- read_xlsx(nri, sheet = "reads_I_II") |>
  clean_names() |>
  mutate(
    sample_id = str_to_lower(str_trim(sample_id)),
    sample_id = if_else(
      sample_id %in% str_to_lower(str_trim(sample_metadata$sample_id)),
      sample_id,
      str_remove(sample_id, "_1$")
    )
  ) |>
  transmute(
    sample_id,
    altitude,
    qc_reads = reads_after_quality_filtering,
    nem_reads = nematode_reads_pident_95,
    prop_nem = nem_reads / qc_reads
  )

sample_metadata <- sample_metadata |>
  left_join(
    reads |> select(sample_id, altitude, qc_reads),
    by = "sample_id"
  )

write_xlsx(sample_metadata, "nematode_sample_metadata.xlsx")
```
Get nematode-only annotations:
```R
#Keep only Nematoda ASVs, keep only needed columns, give both markers the same rank column names (Kingdom ... Species) so they can be compared.

nem_28S <- asv_28S_long |>
  dplyr::filter(tbl_tax_2_phyla == "Nematoda") |>
  dplyr::select(sample_id, asv, count, tbl_tax_2_domain_kingdom, tbl_tax_2_phyla,
                tbl_tax_2_class, tbl_tax_2_order, tbl_tax_2_family,
                tbl_tax_2_genus, tbl_tax_2_species) |>
  dplyr::rename(Kingdom = tbl_tax_2_domain_kingdom, Phylum = tbl_tax_2_phyla,
                Class = tbl_tax_2_class, Order = tbl_tax_2_order,
                Family = tbl_tax_2_family, Genus = tbl_tax_2_genus,
                Species = tbl_tax_2_species)

nem_18S <- asv_18S_long |>
  dplyr::filter(table8_phylum == "Nematoda") |>
  dplyr::select(sample_id, asv, count, table8_domain_kingdom, table8_phylum,
                table8_class, table8_order, table8_family, table8_genus,
                table8_species) |>
  dplyr::rename(Kingdom = table8_domain_kingdom, Phylum = table8_phylum,
                Class = table8_class, Order = table8_order, Family = table8_family,
                Genus = table8_genus, Species = table8_species)

#Tidy the taxonomy strings: trim whitespace, and turn "" and "NA" text into real NA. Otherwise "NA" is counted as a taxon and ends up in shared-taxa lists.
nem_18S <- nem_18S |>
  mutate(across(all_of(c("Class", "Order", "Family", "Genus", "Species")), \(x) x |> str_squish() |> na_if("") |> na_if("NA")))
nem_28S <- nem_28S |>
  mutate(across(all_of(c("Class", "Order", "Family", "Genus", "Species")), \(x) x |> str_squish() |> na_if("") |> na_if("NA")))
```

### Compare 18S and 28S markers

Abundance concordance: Find the % of total ASV count that are from the nematoda phylum:
```R
#For each sample: total reads, reads assigned to phylum Nematoda, and the percentage.

#28S summary
summary_28S <- asv_28S_long |>
  group_by(sample_id) |>
  summarise(
    total_asv_28S = sum(count, na.rm = TRUE),
    nematode_count_28S = sum(count[tbl_tax_2_phyla == "Nematoda"], na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    nematode_percent_28S = if_else(total_asv_28S > 0,
                                   100 * nematode_count_28S / total_asv_28S,
                                   NA_real_)
  )

#18S summary
summary_18S <- asv_18S_long |>
  group_by(sample_id) |>
  summarise(
    total_asv_18S = sum(count, na.rm = TRUE),
    nematode_count_18S = sum(count[table8_phylum == "Nematoda"], na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    nematode_percent_18S = if_else(total_asv_18S > 0,
                                   100 * nematode_count_18S / total_asv_18S,
                                   NA_real_)
  )

sample_metadata <- sample_metadata |>
  select(-any_of(c(names(summary_18S)[-1], names(summary_28S)[-1]))) |>
  left_join(summary_18S, by = "sample_id") |>
  left_join(summary_28S, by = "sample_id")

write_xlsx(sample_metadata, "nematode_sample_metadata.xlsx")
```
Define the paired sample set
```R
MIN_READS <- 0     # min nematode-assigned reads per marker per sample

#Get one row per sample sequenced with BOTH markers, with total reads, nematode reads, and nematode % for each marker. Samples with zero nematode are retained, with 0 (not NA).
depth <- inner_join(
  summary_18S |> select(sample_id, total_18S = total_asv_18S),
  summary_28S |> select(sample_id, total_28S = total_asv_28S),
  by = "sample_id"
) |>
  filter(total_18S > 0, total_28S > 0)

nem_reads_18S <- nem_18S |> group_by(sample_id) |>
  summarise(nem_18S = sum(count, na.rm = TRUE), .groups = "drop")
nem_reads_28S <- nem_28S |> group_by(sample_id) |>
  summarise(nem_28S = sum(count, na.rm = TRUE), .groups = "drop")

depth <- depth |>
  left_join(nem_reads_18S, by = "sample_id") |>
  left_join(nem_reads_28S, by = "sample_id") |>
  mutate(across(c(nem_18S, nem_28S), \(x) replace_na(x, 0)),   # no nematodes = 0 reads
         pct_18S = 100 * nem_18S / total_18S,
         pct_28S = 100 * nem_28S / total_28S) |>
  left_join(sample_metadata |> select(sample_id, location, group, yr), by = "sample_id")

#ids = samples that pass the depth gate in BOTH markers. This is the main analysis set.
ids <- depth |>
  filter(nem_18S >= MIN_READS, nem_28S >= MIN_READS) |>
  pull(sample_id)

cat("Samples analysed:", length(ids), "of", nrow(depth), "paired\n")
print(depth |> filter(sample_id %in% ids) |> count(group))
```
Classification resolution for the two markers for nematode taxa
```R
RANKS     <- c("Class", "Order", "Family", "Genus", "Species")

resolution <- bind_rows(nem_18S |> mutate(marker = "18S"),
                        nem_28S |> mutate(marker = "28S")) |>
  filter(sample_id %in% ids) |>
  pivot_longer(all_of(RANKS), names_to = "Rank", values_to = "taxon") |>
  group_by(marker, Rank) |>
  summarise(pct_reads_classified = 100 * sum(count[!is.na(taxon)]) / sum(count),
            .groups = "drop") |>
  mutate(Rank = factor(Rank, RANKS)) |>
  arrange(marker, Rank)
resolution
#   marker Rank    pct_reads_classified
# 1 18S    Class                  99.5 
# 2 18S    Order                  98.5 
# 3 18S    Family                 96.3 
# 4 18S    Genus                  95.3 
# 5 18S    Species                 3.37
# 6 28S    Class                 100.0 
# 7 28S    Order                  96.1 
# 8 28S    Family                 92.3 
# 9 28S    Genus                  91.0 
#10 28S    Species                77.7 
```

Abundance concordance: are samples / taxa abundant in 18S also abundant in 28S?
```R
#For each filter setting (mr = MIN_REL value) and each rank (rk):
#a. sum reads per sample x taxon, separately for 18S and 28S
#b. convert to within-sample relative abundance (among classified reads)
#c. drop taxa below mr, then renormalise so each sample sums to 1
#d. join the markers into one long table, with a row for EVERY paired sample x EVERY taxon seen at that rank (zeros filled in)
#e. flag samples where both markers resolve >= 1 taxon (`has`)
#f. compute per-sample overlap + Bray-Curtis similarity

MIN_REL   <- 0    # min within-sample relative abundance of a taxon (0 = off)
MIN_REL_GRID   <- sort(unique(c(0, 0.01, 0.05, MIN_REL)))

pair_all       <- list()
per_sample_all <- list()

for (mr in MIN_REL_GRID) {
  key <- as.character(mr)
  pair_all[[key]]       <- list()
  per_sample_all[[key]] <- list()
  for (rk in RANKS) {
    #a-c: 18S relative abundance table (reads without a label at this rank are excluded, so relative abundance is among CLASSIFIED reads)
    t18 <- nem_18S |>
      filter(sample_id %in% depth$sample_id, !is.na(.data[[rk]])) |>
      group_by(sample_id, taxon = .data[[rk]]) |>
      summarise(n_18S = sum(count, na.rm = TRUE), .groups = "drop") |>
      group_by(sample_id) |>
      mutate(rel_18S = n_18S / sum(n_18S)) |>
      filter(rel_18S >= mr) |>
      mutate(rel_18S = rel_18S / sum(rel_18S)) |>    
      ungroup()
    #same for 28S
    t28 <- nem_28S |>
      filter(sample_id %in% depth$sample_id, !is.na(.data[[rk]])) |>
      group_by(sample_id, taxon = .data[[rk]]) |>
      summarise(n_28S = sum(count, na.rm = TRUE), .groups = "drop") |>
      group_by(sample_id) |>
      mutate(rel_28S = n_28S / sum(n_28S)) |>
      filter(rel_28S >= mr) |>
      mutate(rel_28S = rel_28S / sum(rel_28S)) |>
      ungroup()
    #d: join markers. full_join keeps taxa seen by only one marker.
    pr <- full_join(t18, t28, by = c("sample_id", "taxon"))
    if (nrow(pr) == 0) next                          
    pr <- pr |>
      complete(sample_id = depth$sample_id, taxon) |> 
      mutate(across(c(n_18S, n_28S, rel_18S, rel_28S), \(x) replace_na(x, 0))) |>
      #e: `has` = both markers resolve at least one taxon in this sample
      group_by(sample_id) |>
      mutate(has = any(rel_18S > 0) & any(rel_28S > 0)) |>
      ungroup()  
    pair_all[[key]][[rk]] <- pr
    #f: per-sample metrics.
    #   jaccard            = shared taxa / all taxa found by either marker
    #   contain_18S_in_28S = share of 18S taxa that 28S also found
    #   contain_28S_in_18S = share of 28S taxa that 18S also found
    #   bc_similarity      = 1 - Bray-Curtis dissimilarity (1 = identical
    #                        communities, 0 = nothing in common). Both vectors
    #                        sum to 1, so the denominator is 2.
    per_sample_all[[key]][[rk]] <- pr |>
      group_by(sample_id) |>
      summarise(
        n_taxa_18S = sum(rel_18S > 0),
        n_taxa_28S = sum(rel_28S > 0),
        n_shared   = sum(rel_18S > 0 & rel_28S > 0),
        n_union    = sum(rel_18S > 0 | rel_28S > 0),
        jaccard             = n_shared / n_union,
        contain_18S_in_28S  = n_shared / n_taxa_18S,
        contain_28S_in_18S  = n_shared / n_taxa_28S,
        bc_similarity = 1 - sum(abs(rel_18S - rel_28S)) / sum(rel_18S + rel_28S),
        .groups = "drop"
      ) |>
      mutate(has_both = n_taxa_18S > 0 & n_taxa_28S > 0,
             Rank = rk)
  }
}

key_main <- as.character(MIN_REL)

#Taxa never detected by either marker in these samples are dropped, so shared zeros do not inflate correlations or kappa.
main_pair <- list()
for (rk in names(pair_all[[key_main]])) {
  main_pair[[rk]] <- pair_all[[key_main]][[rk]] |>
    filter(sample_id %in% ids) |>
    group_by(taxon) |>
    filter(any(rel_18S > 0 | rel_28S > 0)) |>
    ungroup()
}

#per-sample metrics for `ids`, only samples where both markers resolved something at that rank.
main_ps <- map(per_sample_all[[key_main]],
               \(x) filter(x, sample_id %in% ids, has_both))
```

```R
#Nematode share of reads, per sample, uses ALL paired samples including zeros. Zeros are real information: "this marker saw no nematodes here".
cor.test(depth$pct_18S, depth$pct_28S, method = "spearman", exact = FALSE)
#data:  depth$pct_18S and depth$pct_28S
#S = 40889, p-value = 7.349e-07
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.5207596

ggplot(depth, aes(x = pct_18S, y = pct_28S, colour = group)) +
  geom_point() +
  geom_abline(linetype = 2) +                       # dashed 1:1 line
  scale_x_sqrt() + scale_y_sqrt() +                 # sqrt spreads out the many small values
  labs(x = "18S reads assigned to Nematoda (%)",
       y = "28S reads assigned to Nematoda (%)",
       title = "Nematode share of reads, per sample (all paired samples)") +
  theme_classic()

#Genus level: is a genus abundant in 18S also abundant in 28S? Uses main_pair$Genus: all samples in `ids`, with a row for every sample x genus. `rel_*` is the genus's share of CLASSIFIED nematode reads.
taxon_concordance <- main_pair$Genus |>
  group_by(taxon) |>
  summarise(
    n_samples    = n(),
    n_28S        = sum(rel_28S > 0),                 # samples where 28S detects it
    n_18S        = sum(rel_18S > 0),
    n_both       = sum(rel_18S > 0 & rel_28S > 0),
    mean_rel_18S = mean(rel_18S),
    mean_rel_28S = mean(rel_28S),
    #kappa = agreement on presence/absence beyond chance. Undefined when a
    #marker detects the taxon in none or all samples (no variance).
    kappa = if (n_28S > 0 && n_28S < n_samples && n_18S > 0 && n_18S < n_samples) {
      irr::kappa2(cbind(as.integer(rel_18S > 0), as.integer(rel_28S > 0)))$value
    } else {
      NA_real_
    },
    #rho = do samples with more of this genus in 18S also have more in 28S?
    rho = suppressWarnings(cor(rel_18S, rel_28S, method = "spearman")),
    .groups = "drop"
  ) |>
  arrange(desc(n_28S))
print(taxon_concordance)
#   taxon               n_samples n_28S n_18S n_both mean_rel_18S mean_rel_28S   kappa      rho
# 1 Micoletzkya                80    68    80     68    0.711        0.706     NA      0.567
# 2 Diplogasteroides           80    49    59     46    0.0676       0.0789     0.552  0.716
# 3 Contortylenchus            80    21    28     21    0.210        0.103      0.796  0.900
# 4 Parasitorhabditis          80    21    41     21    0.00431      0.00373    0.506  0.646
# 5 Cryptaphelenchus           80    10    39     10    0.00464      0.0262     0.261  0.547
# 6 Bursaphelenchus            80     8    27      8    0.00217      0.000218   0.358  0.533
# 7 Bunonema                   80     2    15      2    0.000428     0.0000132  0.200  0.342
# 8 Meloidogyne                80     2     0      0    0            0.000388  NA     NA    
# 9 Cephalenchus               80     1     0      0    0            0.00114   NA     NA    
#10 Devibursaphelenchus        80     1     0      0    0            0.00538   NA     NA    
#11 Malenchus                  80     1     0      0    0            0.000238  NA     NA    
#12 Ditylenchus                80     0     6      0    0.0000837    0         NA     NA    
#13 Parasitylenchus            80     0     3      0    0.000243     0         NA     NA 

ggplot(taxon_concordance,
       aes(x = mean_rel_18S, y = mean_rel_28S, label = taxon)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(nudge_y = 0.015, check_overlap = TRUE, size = 3) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Mean genus relative abundance in 18S",
       y = "Mean genus relative abundance in 28S",
       title = "Mean genus abundance concordance between markers") +
  theme_classic()

#Across all genera: do genera abundant in 18S tend to be abundant in 28S? (descriptive only, there are few genera)
cor.test(taxon_concordance$mean_rel_18S, taxon_concordance$mean_rel_28S,
         method = "spearman", exact = FALSE)
#data:  taxon_concordance$mean_rel_18S and taxon_concordance$mean_rel_28S
#S = 148.73, p-value = 0.03327
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.5913995 

ggplot(main_pair$Genus, aes(x = rel_18S, y = rel_28S)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Genus relative abundance in 18S",
       y = "Genus relative abundance in 28S",
       title = "Pooled genus-level abundance concordance") +
  theme_classic()

#Family level fam_has = only samples where both markers resolved >= 1 family. Without this, a sample where 18S could not classify to family would count as "18S: zero".
fam_has <- main_pair$Family |> filter(has)

#including shared zeros (rewards agreeing on absence too)
cor.test(fam_has$rel_18S, fam_has$rel_28S, method = "spearman", exact = FALSE)
#data:  fam_has$rel_18S and fam_has$rel_28S
#S = 71147750, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.7001837 

#only cells where BOTH markers detect the family: when both see it, do they agree on how much?
shared_cells <- fam_has |> filter(rel_18S > 0, rel_28S > 0)
nrow(shared_cells)
#140
cor.test(shared_cells$rel_18S, shared_cells$rel_28S,
         method = "spearman", exact = FALSE)
#data:  shared_cells$rel_18S and shared_cells$rel_28S
#S = 101763, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.7774755

#per family: does a family's abundance track across samples, and which marker reports more? median_log2_28S_over_18S > 0 means 28S reports higher abundance when both detect it; 0 means equal.
per_family <- fam_has |>
  group_by(Family = taxon) |>
  summarise(
    n_either = sum(rel_18S > 0 | rel_28S > 0),
    n_both   = sum(rel_18S > 0 & rel_28S > 0),
    rho = if (n_either >= 5 && sd(rel_18S) > 0 && sd(rel_28S) > 0)
      cor(rel_18S, rel_28S, method = "spearman") else NA_real_,
    median_log2_28S_over_18S = median(log2(rel_28S[rel_18S > 0 & rel_28S > 0] /
                                           rel_18S[rel_18S > 0 & rel_28S > 0])),
    .groups = "drop"
  ) |>
  arrange(desc(n_either))
per_family
#   Family           n_either n_both     rho median_log2_28S_over_18S
# 1 Diplogastridae          75     71  0.646                   0.0186
# 2 Aphelenchoididae        48     19  0.576                  -1.53  
# 3 Rhabditidae             46     27  0.454                   1.18  
# 4 Diplogasteroidea        38      0 NA                      NA     
# 5 Allantonematidae        25     21  0.933                  -1.70  
# 6 Bunonematidae           15      2  0.340                  -1.84  
# 7 Anguinidae               6      0 NA                      NA     
# 8 Enoploidea               5      0 NA                      NA     
# 9 Tylenchidae              5      0 NA                      NA     
#10 Parasitylenchidae        3      0 NA                      NA     
#11 Brevibuccidae            2      0 NA                      NA     
#12 Meloidogynidae           2      0 NA                      NA     
#13 Sphaerularioidea         2      0 NA                      NA     
#14 Aphelenchidae            1      0 NA                      NA     
#15 Pratylenchidae           1      0 NA                      NA  

ggplot(shared_cells, aes(x = rel_18S, y = rel_28S, colour = taxon)) +
  geom_point(alpha = 0.6) +
  geom_abline(linetype = 2) +
  scale_x_sqrt(labels = scales::percent) +
  scale_y_sqrt(labels = scales::percent) +
  labs(x = "18S relative abundance (within nematodes)",
       y = "28S relative abundance (within nematodes)",
       colour = "Family",
       title = "Shared detections, family level") +
  theme_classic()

#abundance relative to ALL reads = (nematode share) x (within-nematode share). If this correlation is weaker than (b), the difference comes from host vs nematode amplification differences between markers, not from community composition.
cells_all <- shared_cells |>
  left_join(depth |> select(sample_id, pct_18S, pct_28S), by = "sample_id") |>
  mutate(all_18S = rel_18S * pct_18S / 100,
         all_28S = rel_28S * pct_28S / 100)
cor.test(cells_all$all_18S, cells_all$all_28S, method = "spearman", exact = FALSE)
#data:  cells_all$all_18S and cells_all$all_28S
#S = 129902, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.7159432
```
Taxonomic concordance: do the two markers identify the same nematode taxa?
```R
#Pooled across the analysed samples, one row per rank.
tc_list <- list()
for (rk in names(main_pair)) {
  d   <- main_pair[[rk]]
  det <- d |>
    group_by(taxon) |>
    summarise(d18 = any(rel_18S > 0), d28 = any(rel_28S > 0), .groups = "drop")
  dh  <- d |> filter(has)
  
  tc_list[[rk]] <- tibble(
    Rank       = rk,
    n_samples  = n_distinct(dh$sample_id),
    n_taxa_18S = sum(det$d18),
    n_taxa_28S = sum(det$d28),
    n_shared   = sum(det$d18 & det$d28),
    n_union    = sum(det$d18 | det$d28),
    jaccard    = n_shared / n_union,
    prop_18S_in_28S = n_shared / n_taxa_18S,
    prop_28S_in_18S = n_shared / n_taxa_28S,
    kappa_presence  = if (n_distinct(dh$sample_id) > 1)
      kappa2(cbind(as.integer(dh$rel_18S > 0), as.integer(dh$rel_28S > 0)))$value
    else NA_real_
  )
}
taxonomic_concordance <- bind_rows(tc_list)
taxonomic_concordance
#  Rank   n_samples n_taxa_18S n_taxa_28S n_shared n_union jaccard prop_18S_in_28S prop_28S_in_18S
#1 Class         75          2          2        1       3   0.333           0.5             0.5  
#2 Order         75          4          5        3       6   0.5             0.75            0.6  
#3 Family        75         11          9        5      15   0.333           0.455           0.556
#4 Genus         74          9         11        7      13   0.538           0.778           0.636
#5 Speci…        60          5          5        1       9   0.111           0.2             0.2

#Within each sample.
per_sample_concordance <- bind_rows(main_ps)
per_sample_summary <- per_sample_concordance |>
  mutate(Rank = factor(Rank, RANKS)) |>
  group_by(Rank) |>
  summarise(n = n(),
            median_jaccard = median(jaccard),
            median_contain_18S_in_28S = median(contain_18S_in_28S),
            .groups = "drop")
per_sample_summary
#  Rank        n median_jaccard median_contain_18S_in_28S
#1 Class      75          1                         1    
#2 Order      75          0.6                       0.667
#3 Family     75           0.5                      0.5   
#4 Genus      74          0.6                       0.633
#5 Species    60          0.25                      0.5 

#Genus detection frequency by marker (number of samples each marker detects it in)
ggplot(taxon_concordance, aes(x = reorder(taxon, n_28S), y = n_28S)) +
  geom_col(aes(fill = "28S")) +
  geom_col(aes(y = n_18S, fill = "18S"), alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("18S" = "steelblue", "28S" = "darkorange"),
                    name = "Marker") +
  labs(x = "Genus", y = "Number of samples detected",
       title = "Genus detection frequency by marker") +
  theme_classic()

#Which genera are shared / marker-specific?
genus_detect <- main_pair$Genus |>
  group_by(taxon) |>
  summarise(d18 = any(rel_18S > 0), d28 = any(rel_28S > 0), .groups = "drop")
list(shared   = genus_detect$taxon[genus_detect$d18 & genus_detect$d28],
     only_18S = genus_detect$taxon[genus_detect$d18 & !genus_detect$d28],
     only_28S = genus_detect$taxon[!genus_detect$d18 & genus_detect$d28])

#shared
# "Bunonema"          "Bursaphelenchus"   "Contortylenchus"   "Cryptaphelenchus" 
# "Diplogasteroides"  "Micoletzkya"       "Parasitorhabditis"

#only_18S
# "Ditylenchus"     "Parasitylenchus"

#only_28S
# "Cephalenchus"        "Devibursaphelenchus" "Malenchus"           "Meloidogyne"  
```
Detection sensitivity: how often does 18S see what 28S sees?
```R
#Take every sample x genus cell that 28S detects and ask whether 18S also detects it. Rhis includes samples where 18S could not classify to genus at all, which is a real limit of 18S rather than a missing sample.
detection <- main_pair$Genus |>
  filter(rel_28S > 0) |>
  mutate(det_18S = rel_18S > 0)

#By genus
detection_by_taxon <- detection |>
  group_by(taxon) |>
  summarise(n_28S = n(),
            prop_detected_18S = mean(det_18S),
            median_rel_28S = median(rel_28S),
            .groups = "drop") |>
  arrange(desc(n_28S))
detection_by_taxon
#   taxon               n_28S prop_detected_18S median_rel_28S
# 1 Micoletzkya            68             1           0.962   
# 2 Diplogasteroides       49             0.939       0.0126  
# 3 Contortylenchus        21             1           0.299   
# 4 Parasitorhabditis      21             1           0.00111 
# 5 Cryptaphelenchus       10             1           0.00339 
# 6 Bursaphelenchus         8             1           0.000787
# 7 Bunonema                2             1           0.000527
# 8 Meloidogyne             2             0           0.0155  
# 9 Cephalenchus            1             0           0.0912  
#10 Devibursaphelenchus     1             0           0.431   
#11 Malenchus               1             0           0.0190 

#Bin how abundant the genus is in 28S. If 18S misses the genera 28S finds most abundant, the gap is more notable.
detection_by_bin <- detection |>
  mutate(bin = cut(rel_28S, c(0, .01, .05, .25, 1), include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(n = n(), prop_detected_18S = mean(det_18S), .groups = "drop")
detection_by_bin
#  bin             n prop_detected_18S
#1 [0,0.01]       60             0.967
#2 (0.01,0.05]    17             0.824
#3 (0.05,0.25]    20             0.95 
#4 (0.25,1]       87             0.977
```
Sample-level concordance: do the markers recover similar communities?
```R
#Bray-Curtis similarity per sample at genus and family level 
similarity_summary <- bind_rows(main_ps$Genus, main_ps$Family) |>
  group_by(Rank) |>
  summarise(n_samples = n(),
            median = median(bc_similarity, na.rm = TRUE),
            mean   = mean(bc_similarity, na.rm = TRUE),
            IQR    = IQR(bc_similarity, na.rm = TRUE),
            .groups = "drop")
print(similarity_summary)
#  Rank   n_samples  median  mean    IQR
#1 Family        75  0.947 0.791 0.411
#2 Genus         74 0.897  0.772 0.401

#Family-level distribution 
fam_sim <- main_ps$Family |>
  left_join(depth |> select(sample_id, group, nem_18S, nem_28S), by = "sample_id")
summary(fam_sim$bc_similarity)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000718 0.579294 0.946732 0.791477 0.989931 1.000000   

ggplot(bind_rows(main_ps$Genus, main_ps$Family), aes(x = bc_similarity)) +
  geom_histogram(bins = 20) +
  facet_wrap(~Rank) +
  labs(x = "Bray-Curtis similarity between markers",
       y = "Number of samples",
       title = "18S-28S sample-level community concordance") +
  theme_classic()

#depth effect?
ggplot(fam_sim, aes(x = pmin(nem_18S, nem_28S), y = bc_similarity, colour = group)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Nematode reads (lower of the two markers)", y = "BC similarity (family)") +
  theme_classic()

#Mantel test: do the markers agree on which samples resemble each other? Builds a sample x taxon matrix per marker, turns each into a Bray-Curtis distance matrix between samples, and correlates the two distance matrices. A significant positive r means among-sample community structure is shared.
mantel_list <- list()
for (rk in c("Family", "Genus")) {
  d <- main_pair[[rk]] |> filter(has)
  if (n_distinct(d$sample_id) < 5) next                
  M18 <- d |> select(sample_id, taxon, rel_18S) |>
    pivot_wider(names_from = taxon, values_from = rel_18S) |>
    column_to_rownames("sample_id") |> as.matrix()
  M28 <- d |> select(sample_id, taxon, rel_28S) |>
    pivot_wider(names_from = taxon, values_from = rel_28S) |>
    column_to_rownames("sample_id") |> as.matrix()
  m <- vegan::mantel(vegdist(M18, "bray"), vegdist(M28, "bray"),
                     method = "spearman", permutations = 9999)
  mantel_list[[rk]] <- tibble(Rank = rk, n = nrow(M18),
                              mantel_r = m$statistic, mantel_p = m$signif)
}
mantel_tbl <- bind_rows(mantel_list)
mantel_tbl
#  Rank       n mantel_r mantel_p
#1 Family    75    0.606   0.0001
#2 Genus     74    0.556   0.0001

#Similarity versus read depth
depth_effect <- bind_rows(per_sample_all[[key_main]]$Genus,
                          per_sample_all[[key_main]]$Family) |>
  filter(!is.na(bc_similarity)) |>
  left_join(depth |> select(sample_id, group, nem_18S, nem_28S), by = "sample_id") |>
  mutate(min_reads = pmin(nem_18S, nem_28S))

depth_effect |>
  group_by(Rank) |>
  summarise(n = n(),
            rho = cor(min_reads, bc_similarity, method = "spearman"),
            .groups = "drop")
#  Rank       n   rho
#1 Family    80 0.0272
#2 Genus     80 0.0426

ggplot(depth_effect, aes(x = min_reads, y = bc_similarity, colour = group)) +
  geom_point() +
  geom_vline(xintercept = MIN_READS, linetype = 2) +    # the depth gate
  scale_x_log10() +
  facet_wrap(~Rank) +
  labs(x = "Nematode reads (lower of the two markers)",
       y = "Bray-Curtis similarity between markers",
       colour = "Group") +
  theme_classic()
```
Sensitivity to low depth/abundance taxa
```R
MIN_READS_GRID <- sort(unique(c(0, 50, 100, 500, MIN_READS)))
MIN_REL_GRID   <- sort(unique(c(0, 0.01, 0.05, MIN_REL)))

sens_grid <- expand_grid(min_reads = MIN_READS_GRID, min_rel = MIN_REL_GRID)
sens_rows <- list()

for (i in seq_len(nrow(sens_grid))) {
  mn  <- sens_grid$min_reads[i]
  mr  <- sens_grid$min_rel[i]
  key <- as.character(mr)
  
  ids_s <- depth |> filter(nem_18S >= mn, nem_28S >= mn) |> pull(sample_id)
  
  # family table for this setting, restricted like main_pair
  fam_s <- pair_all[[key]]$Family |>
    filter(sample_id %in% ids_s) |>
    group_by(taxon) |>
    filter(any(rel_18S > 0 | rel_28S > 0)) |>
    ungroup()
  fam_det <- fam_s |>
    group_by(taxon) |>
    summarise(d18 = any(rel_18S > 0), d28 = any(rel_28S > 0), .groups = "drop")
  fam_h <- fam_s |> filter(has)
  
  ps_fam <- per_sample_all[[key]]$Family |> filter(sample_id %in% ids_s, has_both)
  ps_gen <- per_sample_all[[key]]$Genus  |> filter(sample_id %in% ids_s, has_both)
  
  # Mantel (family) with fewer permutations, since it is run for every setting
  mantel_fam <- NA_real_
  if (n_distinct(fam_h$sample_id) >= 5) {
    M18 <- fam_h |> select(sample_id, taxon, rel_18S) |>
      pivot_wider(names_from = taxon, values_from = rel_18S) |>
      column_to_rownames("sample_id") |> as.matrix()
    M28 <- fam_h |> select(sample_id, taxon, rel_28S) |>
      pivot_wider(names_from = taxon, values_from = rel_28S) |>
      column_to_rownames("sample_id") |> as.matrix()
    mantel_fam <- vegan::mantel(vegdist(M18, "bray"), vegdist(M28, "bray"),
                                method = "spearman", permutations = 999)$statistic
  }
  
  both <- fam_h$rel_18S > 0 & fam_h$rel_28S > 0
  
  sens_rows[[i]] <- tibble(
    min_reads = mn, min_rel = mr,
    n_samples = length(ids_s), n_paired_total = nrow(depth),
    fam_prop_18S_in_28S = sum(fam_det$d18 & fam_det$d28) / sum(fam_det$d18),
    fam_median_bc_sim   = median(ps_fam$bc_similarity),
    gen_median_bc_sim   = median(ps_gen$bc_similarity),
    mantel_family       = mantel_fam,
    rho_all_cells       = cor(fam_h$rel_18S, fam_h$rel_28S, method = "spearman"),
    rho_shared_cells    = cor(fam_h$rel_18S[both], fam_h$rel_28S[both], method = "spearman"),
    n_shared_cells      = sum(both)
  )
}
sens <- bind_rows(sens_rows)
sens
#   min_reads min_rel n_samples n_paired_total fam_prop_18S_in_28S fam_median_bc_sim gen_median_bc_sim mantel_family rho_all_cells rho_shared_cells n_shared_cells
# 1         0    0           80             80               0.455             0.947             0.897         0.606         0.700            0.777            140
# 2         0    0.01        80             80               0.5               0.950             0.899         0.606         0.770            0.626             94
# 3         0    0.05        80             80               0.8               1                 0.898         0.610         0.864            0.641             86
# 4        50    0           64             80               0.455             0.946             0.883         0.702         0.707            0.832            129
# 5        50    0.01        64             80               0.571             0.949             0.887         0.695         0.806            0.640             85
# 6        50    0.05        64             80               0.8               1                 0.887         0.713         0.904            0.615             78
# 7       100    0           60             80               0.455             0.948             0.883         0.736         0.713            0.839            124
# 8       100    0.01        60             80               0.571             0.949             0.887         0.730         0.811            0.655             80
# 9       100    0.05        60             80               0.8               1                 0.887         0.737         0.910            0.608             74
#10       500    0           54             80               0.455             0.955             0.897         0.739         0.721            0.841            116
#11       500    0.01        54             80               0.8               0.952             0.899         0.732         0.831            0.651             73
#12       500    0.05        54             80               0.75              1                 0.898         0.733         0.909            0.588             67 

#Values barely move across rows, therefore they are robust to the filters. 
```
Export
```R
write_xlsx(list(
  depth                  = depth,
  resolution             = resolution,
  taxonomic_concordance  = taxonomic_concordance,
  per_sample_concordance = per_sample_concordance,
  per_sample_summary     = per_sample_summary,
  taxon_concordance      = taxon_concordance,
  per_family             = per_family,
  detection_by_taxon     = detection_by_taxon,
  detection_by_bin       = detection_by_bin,
  similarity_summary     = similarity_summary,
  mantel                 = mantel_tbl,
  sensitivity            = sens
), "marker_concordance_18S_28S.xlsx")
```
### Diversity

#### All
Create phyloseq objects

```R
#28S
meta <- sample_metadata |>
  mutate(sample_key = str_to_lower(str_trim(sample_id)))

df28S <- asv_28S_long |>
  dplyr::select(sample_id, asv, count, tbl_tax_2_domain_kingdom, tbl_tax_2_phyla,
                tbl_tax_2_class, tbl_tax_2_order, tbl_tax_2_family,
                tbl_tax_2_genus, tbl_tax_2_species) |>
  dplyr::rename(Kingdom = tbl_tax_2_domain_kingdom, Phylum = tbl_tax_2_phyla,
                Class = tbl_tax_2_class, Order = tbl_tax_2_order,
                Family = tbl_tax_2_family, Genus = tbl_tax_2_genus,
                Species = tbl_tax_2_species) |>
  mutate(across(c(Class, Order, Family, Genus, Species), ~ .x |> str_squish() |> na_if("") |> na_if("NA")))

otu_wide <- df28S |> dplyr::group_by(sample_id, asv) |> dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop"
  ) |> tidyr::pivot_wider(names_from = asv, values_from = count, values_fill = 0)

otu <- as.matrix(otu_wide[, -1, drop = FALSE])
rownames(otu) <- otu_wide$sample_id
storage.mode(otu) <- "numeric"
otu <- t(otu)

meta_28 <- as.data.frame(meta)
rownames(meta_28) <- meta_28$sample_key

colnames(otu) <- str_to_lower(str_trim(colnames(otu)))

common_samples <- intersect(colnames(otu), rownames(meta_28))

otu <- otu[, common_samples, drop = FALSE]
meta_28 <- meta_28[common_samples, , drop = FALSE]

tax_split <- df28S |> dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species) |> dplyr::distinct()

if (anyDuplicated(tax_split$asv)) {
  stop("Some 28S ASVs have conflicting taxonomy assignments.")
}

tax_split <- as.data.frame(tax_split)
rownames(tax_split) <- tax_split$asv
tax_split$asv <- NULL

tax_split <- tax_split[rownames(otu), , drop = FALSE]

OTU <- otu_table(otu, taxa_are_rows = TRUE)
SAM <- sample_data(meta_28)
TAX <- tax_table(as.matrix(tax_split))
ps_28 <- phyloseq(OTU, SAM, TAX)
 
 
#18S
df18S <- asv_18S_long |>
  dplyr::select(sample_id, asv, count, table8_domain_kingdom, table8_phylum,
                table8_class, table8_order, table8_family, table8_genus,
                table8_species) |>
  dplyr::rename(Kingdom = table8_domain_kingdom, Phylum = table8_phylum,
                Class = table8_class, Order = table8_order, Family = table8_family,
                Genus = table8_genus, Species = table8_species) |>
  mutate(across(c(Class, Order, Family, Genus, Species), ~ .x |> str_squish() |> na_if("") |> na_if("NA")))

otu_wide <- df18S |> dplyr::group_by(sample_id, asv) |> dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop"
  ) |> tidyr::pivot_wider(names_from = asv, values_from = count, values_fill = 0)

otu <- as.matrix(otu_wide[, -1, drop = FALSE])
rownames(otu) <- otu_wide$sample_id
storage.mode(otu) <- "numeric"
otu <- t(otu)

meta_18 <- as.data.frame(meta)
rownames(meta_18) <- meta_18$sample_key

colnames(otu) <- str_to_lower(str_trim(colnames(otu)))

common_samples <- intersect(colnames(otu), rownames(meta_18))

otu <- otu[, common_samples, drop = FALSE]
meta_18 <- meta_18[common_samples, , drop = FALSE]

tax_split <- df18S |> dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species) |> dplyr::distinct()

if (anyDuplicated(tax_split$asv)) {
  stop("Some 18S ASVs have conflicting taxonomy assignments.")
}

tax_split <- as.data.frame(tax_split)
rownames(tax_split) <- tax_split$asv
tax_split$asv <- NULL

tax_split <- tax_split[rownames(otu), , drop = FALSE]

OTU <- otu_table(otu, taxa_are_rows = TRUE)
SAM <- sample_data(meta_18)
TAX <- tax_table(as.matrix(tax_split))
ps_18 <- phyloseq(OTU, SAM, TAX)
```
```R
get_alpha <- function(ps_obj, meta_obj) {
  alpha <- estimate_richness(ps_obj,
                             measures = c("Shannon", "Simpson", "Chao1", "Observed"))
  rownames(alpha) <- gsub("\\.", "-", rownames(alpha))
  alpha$Sample <- rownames(alpha)
  meta_obj$Sample <- rownames(meta_obj)
  df <- merge(alpha, meta_obj, by = "Sample")
  return(df)
}

meta_f <- meta_28
ps_f <- prune_samples(rownames(meta_f), ps_28)

alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = group),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 80, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")

#### innichen only
meta_f <- meta_28[meta_28$location == "Innichen", ]
ps_f <- prune_samples(rownames(meta_f), ps_28)

alpha_df <- get_alpha(ps_f, meta_f)

alpha_df <- alpha_df %>%
  mutate(yr = as.factor(yr))
  
ggplot(alpha_df, aes(x = yr, y = Shannon, fill = yr)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = group),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 80, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")
```
```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = group, shape = region)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Sample),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

ggplot(pcoa_df, aes(PC1, PC2, color = yr, shape = location)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

#### innichen only
pcoa_df <- pcoa_df %>%
  mutate(yr = as.factor(yr))

ggplot(pcoa_df, aes(PC1, PC2, color = group, shape = yr)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )
```
```R
meta_f <- meta_18
ps_f <- prune_samples(rownames(meta_f), ps_18)

alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = group),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 80, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")
```
```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = group, shape = region)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Sample),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )
```
#### Nematoda
Create phyloseq objects
```R
ps_18_nematoda <- subset_taxa(
  ps_18,
  Phylum == "Nematoda"
)

ps_28_nematoda <- subset_taxa(
  ps_28,
  Phylum == "Nematoda"
)
```
```R
get_alpha <- function(ps_obj, meta_obj) {
  alpha <- estimate_richness(ps_obj,
                             measures = c("Shannon", "Simpson", "Chao1", "Observed"))
  rownames(alpha) <- gsub("\\.", "-", rownames(alpha))
  alpha$Sample <- rownames(alpha)
  meta_obj$Sample <- rownames(meta_obj)
  df <- merge(alpha, meta_obj, by = "Sample")
  return(df)
}

meta_f <- meta_28
ps_f <- prune_samples(rownames(meta_f), ps_28_nematoda)
ps_f <- prune_samples(sample_sums(ps_f) > 0, ps_f)

alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = group),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 80, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")

#### innichen only
meta_f <- meta_28[meta_28$location == "Innichen", ]
ps_f <- prune_samples(rownames(meta_f), ps_28_nematoda)
ps_f <- prune_samples(sample_sums(ps_f) > 0, ps_f)

alpha_df <- get_alpha(ps_f, meta_f)

alpha_df <- alpha_df %>%
  mutate(yr = as.factor(yr))
  
ggplot(alpha_df, aes(x = yr, y = Shannon, fill = yr)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = group),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 80, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")
```
```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = group, shape = region)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Sample),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

ggplot(pcoa_df, aes(PC1, PC2, color = group, shape = region)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

#### innichen only
pcoa_df <- pcoa_df %>%
  mutate(yr = as.factor(yr))

ggplot(pcoa_df, aes(PC1, PC2, color = group, shape = yr)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )
```
```R
meta_f <- meta_18
ps_f <- prune_samples(rownames(meta_f), ps_18_nematoda)
ps_f <- prune_samples(sample_sums(ps_f) > 0, ps_f)

alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = group),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 80, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")
```
```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = group, shape = region)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Sample),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )
```







































Create phyloseq objects
```R
#28S
otu_28S <- nem_28S |>
  group_by(sample_id, asv) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    names_from = asv,
    values_from = count,
    values_fill = 0
  )

otu_mat_28S <- as.matrix(otu_28S[, -1, drop = FALSE])
rownames(otu_mat_28S) <- otu_28S$sample_id
storage.mode(otu_mat_28S) <- "numeric"

tax_28S <- nem_28S |>
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species) |>
  distinct()

if (anyDuplicated(tax_28S$asv)) {
  stop("Some 28S ASVs have conflicting taxonomy assignments.")
}

tax_mat_28S <- as.matrix(tax_28S[, -1, drop = FALSE])
rownames(tax_mat_28S) <- tax_28S$asv

ps_28S <- phyloseq(
  otu_table(otu_mat_28S, taxa_are_rows = FALSE),
  tax_table(tax_mat_28S)
)


#18S
otu_18S <- nem_18S |>
  group_by(sample_id, asv) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    names_from = asv,
    values_from = count,
    values_fill = 0
  )

otu_mat_18S <- as.matrix(otu_18S[, -1, drop = FALSE])
rownames(otu_mat_18S) <- otu_18S$sample_id
storage.mode(otu_mat_18S) <- "numeric"

tax_18S <- nem_18S |>
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species) |>
  distinct()

if (anyDuplicated(tax_18S$asv)) {
  stop("Some 18S ASVs have conflicting taxonomy assignments.")
}

tax_mat_18S <- as.matrix(tax_18S[, -1, drop = FALSE])
rownames(tax_mat_18S) <- tax_18S$asv

ps_18S <- phyloseq(
  otu_table(otu_mat_18S, taxa_are_rows = FALSE),
  tax_table(tax_mat_18S)
)
```


1. Pairwise presence/absence models with the event as a random effect. This is the most transparent test, for a few biologically motivated pairs (phoretic × phoretic, parasite × parasite, M. buetschlii × each parasite):
```r
   library(glmmTMB)
   # pa: beetle x taxon 0/1 (presence = reads >= your minimum, e.g. >=10 reads AND >=0.1%)
   # d: one row per beetle with sample_id, group (event), nem_reads, and pa columns
   fit_pair <- function(y, x) {
     f <- as.formula(paste(y, "~", x, "+ log(nem_reads) + (1|group)"))
     m <- glmmTMB(f, data = d, family = binomial)
     broom.mixed::tidy(m, effects = "fixed") |> filter(term == x)
   }
   pairs <- tribble(~y, ~x,
     "Contortylenchus", "Micoletzkya_buetschlii",
     "Parasitorhabditis_obtusa", "Micoletzkya_buetschlii",
     "Contortylenchus", "Parasitorhabditis_obtusa",
     "Diplogasteroides_halleri", "Micoletzkya_buetschlii")
   pmap_dfr(pairs, fit_pair, .id = "pair") |> mutate(p_adj = p.adjust(p.value, "holm"))
```

A significant negative coefficient means the taxa segregate beyond what event and depth explain. Compositional dependence is much weaker with presence/absence, though depth-dependent detection remains, hence the depth term.

2. A joint species distribution model (HMSC or ecoCopula). These estimate a residual association matrix between species after accounting for covariates (event, depth, altitude) and are the closest equivalent to your network, built for ~10 species and ~200 sites. For HMSC, use presence/absence with a probit model, event as a random effect, and log(nem_reads) as a covariate.
3. A stratified null-model test. Run vegan::oecosimu (C-score, curveball) separately within each sampling event, or permute within events, to test whether the overall pattern is more segregated than chance.


```R
# drop list from the rule, not by hand
drop_rule <- reads |> filter(prop_nem < 1e-4) |> pull(sample_id)
stopifnot(length(drop_rule) == 18)          # report says 18
setdiff(drop_rule, read_xlsx(veg, "samples_to_drop")$sample_id)  # -> "RV141"

meta <- meta |> left_join(reads, by = "sample_id") |>
  mutate(keep = !sample_id %in% drop_rule)

blast_long <- read_xlsx(nri, "28S") |>
  mutate(tax = clean_tax(species), family = str_squish(family)) |>
  select(sample_id, tax, counts, total_nematode_reads)

blast_counts <- blast_long |> select(sample_id, tax, counts) |>
  pivot_wider(names_from = tax, values_from = counts, values_fill = 0) |>
  column_to_rownames("sample_id")

# Full multi-classifier comparison (wide), if you want to re-derive concordance stats:
merged_28S_wide  <- read_xlsx(path, sheet = "merged_ASV_classifications_28S",
                               skip = 1, na = c("", "NA")) |> clean_names()
```
Build a clean, analysis-ready table
```R
# Long table: one row per Sample x ASV, with manual taxonomy attached, nematodes only
nem_28S <- asv_28S_long |>
  filter(domain_kingdom == "metazoa", phyla == "Nematoda") |>   # adjust col names after clean_names()
  select(sample, asv, count, genus, species, justification)

# Parse sample names into Location / Year / Month using the same logic as your BLAST pipeline
# e.g. "28S_D2A_RV065" -> alpine run; "28S_Ab1_1" -> Abetone
# Best practice: build an explicit sample metadata table (one row per sample) with
# Location, Country, Year, Month, Altitude — matching Table 1 of the PDF — and join it in,
# rather than parsing strings inline every time.
sample_metadata <- read_xlsx("sample_metadata.xlsx")   # you'll need to create/export this once

nem_28S <- nem_28S |> left_join(sample_metadata, by = "sample")
```
From long counts to a community matrix (for vegan/PERMANOVA/NMDS)

Concordance checks (ASV pipeline vs BLAST/BASTA report pipeline)

```R
ids <- meta |> filter(keep) |> pull(sample_id)
align <- function(m, ids, cols) { m <- m[ids, , drop = FALSE]
  miss <- setdiff(cols, colnames(m)); m[miss] <- 0; m[, cols, drop = FALSE] }
cols <- union(colnames(asv_genus), colnames(blast_genus))
A <- align(asv_genus,   ids, cols) |> relab(thr = 0.01)   # same filter on both
B <- align(blast_genus, ids, cols) |> relab(thr = 0.01)

# presence/absence agreement and abundance correlation, per core genus
core <- c("Micoletzkya","Diplogasteroides","Contortylenchus","Parasitorhabditis")
map_dfr(core, \(g) tibble(genus = g,
  kappa = kappa2(cbind(A[, g] > 0, B[, g] > 0))$value,
  rho   = cor(A[, g], B[, g], method = "spearman", use = "complete.obs")))

# community-level
ok <- complete.cases(A, B)
dA <- vegdist(A[ok, ], "bray"); dB <- vegdist(B[ok, ], "bray")
mantel(dA, dB, permutations = 9999)
protest(metaMDS(A[ok, ], "bray", k = 2, trymax = 200),
        metaMDS(B[ok, ], "bray", k = 2, trymax = 200), permutations = 9999)

# same model on both pipelines, side by side
md <- meta |> filter(sample_id %in% ids[ok])
perm <- \(d) adonis2(d ~ group, data = md, permutations = 9999, by = "margin")
list(ASV = perm(dA), BLAST = perm(dB))
````



























## Other 

Classify ASVs with Kraken2 ?

Kraken2/Bracken with a SILVA-based database. Kraken2 uses exact k-mer matches against a database - closest to BOLD. Taxonomic precision depends strongly on database quality, works better with large databases of thousands of reads/ASVs.

VSEARCH/BLAST + SILVA fasta. VSEARCH / BLAST against SILVA/PR2 (very common)
