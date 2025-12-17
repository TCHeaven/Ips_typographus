# Reanalysis of V. Rau <i>Ips typographus</i> / nematode 18S and 28S metabarcoding dataset

Veronica Rau originally generated a dataset of 18S and 28S metabarcoding reads as part of a project to investigate, primarily, the diversity of nematodes associated with <i>Ipy typographus</i> in the S. Tyrol region using molecular methods, as previous identification has been done via morphology only.

D2A and D3B primers were used to target the 28S region: ACAAGTACCGTGAGGGAAAGTTG + TCGGAAGGAACCAGCTACTA
NEM (forward and reverse) primers were used to target the 18S region: GCAAGTCTGGTGCCAGCAGC + CCGTGTTGAGTCAAATTAAG

These primers were chosen as they target the nematode 18 and 28S regions. However, excluding reads aligning to the host <i>I. typographus</i> genome, >25% of reads were classified to nematoda via BLASTP search, and <50% in any given sample. This left many reads unassigned.

We are interested to know what other species are represented in the dataset amongst the non-nematoda, non-Ips reads. Maja Fluch investigated further by using the BOLD (Barcode Of Life Data Systems) database to assign identity to the reads as an alternative approach to BLASTP against the NCBI nt database (Veronica's method). However, the results of this analysis remain confusing as many taxa are assigned that seem biologically implausible.

My aim is the draw out some clearer patterns from the data - if possible. Below are documented analysis and commands used with this dataset. All commands were executed on the ScientificNet HPC Cluster accessed at the University of Bozen-Bolzano from the directory /data/users/theaven/nematode_project, unless stated otherwise.

## Contents

1. [Collecting Data](#2)
2. [Quality Control](#3)
	2.1 [FastQC](#7)
	2.2 [Cutadapt](#8)
3. [Placeholder](#4)
4. [Placeholder](#5)
5. [Placeholder](#6)

DADA2 models and removes error patterns from Illumina data. Deblur, unoise/usearch/vsearch are alternative denoising tools.
Trim off ends where quality drops (e.g. truncate at Q < 20–25). Discard reads with; 1 expected error (EE) or a max EE threshold (e.g. 1–2 per read), Any ambiguous base (N). Most denoisers have built-in recommended parameters (e.g. maxEE in DADA2). DADA can be used inside of QIIME2 - In this case, DADA2 is already your ASV engine, and QIIME2 is just doing taxonomy and stats on the ASVs. Or you can run DADA2 and then export the ASVs to QIIME2.

Merge forward and reverse reads in the overlapping region. When F and R disagree, the merge algorithm uses quality scores to pick the most likely base or discards the read.

Remove chimeras: Use chimera detection (uchime, DADA2’s removeBimeraDenovo, VSEARCH, etc.) Chimeras often: Are low abundance. Are composed of 5’ from one abundant ASV and 3’ from another. Filtering them reduces spurious diversity.

Because systematic errors often appear at low counts: Drop ASVs that: Occur in only one sample, and Are below a minimal count or relative abundance threshold (e.g. <10 reads or <0.1% of sample reads). Or require ASVs to be present in ≥2 samples to be retained. This is a blunt tool (you can throw away real rare things), but in practice it helps a lot for error pruning.

1.Denoise with DADA2 or Deblur. Run a denoiser (DADA2, Deblur, UNOISE) to get ASVs. Reads will usually collapse into far fewer ASVs.

2.Remove chimeras

3.Apply minimum-abundance / prevalence filters

4.Classify ASVs with Kraken2 (or BLAST)

5.Remove ASVs classified outside expected taxa

Kraken2/Bracken with a SILVA-based database. Kraken2 uses exact k-mer matches against a database - closest to BOLD. Taxonomic precision depends strongly on database quality, works better with large databases of thousands of reads/ASVs.

SILVA SSU and LSU - BLAST locally
SILVA ACT - SILVA’s ACT service allows you to upload your own sequences and have them aligned, taxonomically classified, and placed into a reference phylogenetic tree. This is the official and recommended way to classify your 18S/28S sequences using SILVA. You do not need to submit sequences to the SILVA database itself; you simply upload them to the ACT tool for annotation. For large numbers of sequences, ACT can be slow; a local BLAST or alignment against SILVA may be faster. SILVA releases can be downloaded from the SILVA site (including pre-formatted QIIME2/DADA2 training sets

BOLD - Boldigger - k-mer / distance-based matching - its own alignment and similarity scoring (optimized for COI) - BIN system (Barcode Index Number clusters—computed with algorithms not based on BLAST)

If you need high-throughput classification, many researchers use:
QIIME2 + SILVA classifier. QIIME2 is a naive Bayes classifier - Learns k-mer frequency patterns. Handles short, conserved sequences better than BLAST and is better for assigning ambiguous reads, handling subtle variations in conserved regions, and avoiding overconfident species matches. Gives probability-based taxonomic assignments. Train a classifier on SILVA/PR2 trimmed to your primer region. Run qiime feature-classifier classify-sklearn. Get per-ASV taxonomy (up to genus/species, if present in the DB). QIIME2 Naïve Bayes needs a database trimmed to your exact primer region.

VSEARCH/BLAST + SILVA fasta. VSEARCH / BLAST against SILVA/PR2 (very common)

SINA aligner (SILVA’s standalone alignment tool)

IDTAXA (DECIPHER package in R) - probabilistic sequence classification, model-based learning - usually more accurate and conservative than BLAST or naive Bayes. IDTAXA does not require trimming and handles full-length sequences correctly. IDTAXA is slower than QIIME2 NB, but usually much more accurate for 18S/28S metabarcoding.

BLAST is alignment based, limited by 'best hit' interpretation, and prone to misidentifying short or conserved seqeunces. For eDNA/metabarcoding, BLAST is often too literal — it finds the closest sequence, even if it’s wrong.

NemaBase - “nematode-optimized” subset of 18S rRNA reference sequences. Source: all nematode 18S sequences from SILVA (v111 and v138). Curation: cleaning up taxonomy (standardizing taxonomic ranks, updating according to accepted nematode classification, removing suspicious or erroneous entries) to produce a cleaner, well-annotated dataset. Best for soil, freshwater or terrestrial nematode communities where sequenced nematodes are common in SILVA and common metabarcoding primers are used.

Charrier et al. 2024 rRNA-cistron 18S DB - Four databases; 18S, 28S, ITS-1, and ITS-1-5.8S-ITS-2. Curation: Built with the markerDB pipeline (systematic mining and filtering of GenBank) and deduplicated to full-length representatives. Has maximum species coverage across Nematoda. Best for wildlife parasites, marine/free-living nematodes, or very broad community surveys where phylum-wide coverage matters.

NemaTaxa - Built explicitly from sequences that match the NF1 / 18Sr2b primer pair (the standard soil nematode metabarcoding marker). Source: Nematode 18S sequences pulled from NCBI + SILVA v132. Curation: Manual trimming of taxonomy to classical Linnaean ranks (kingdom → genus). Missing ranks filled manually; inconsistent strings cleaned. Good for the the NF1 / 18Sr2b primer pair soil/agricultural metabarcoding. Very clean taxonomy strings suitable for QIIME2.

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

## Quality Control <a name="3"></a>

#### FastQC  <a name="7"></a>
The raw sequence reads were subjected to a quality control check using FastQC.
```bash
screen -S nematode
module load anaconda3
for ReadDir in $(ls -d /data/users/theaven/nematode_project/raw_data/*S/*); do
	Task=FastQC
	ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
	OutDir="$(dirname "$ReadDir")/"$Task""
	ExpectedOutput="$OutDir"/$(basename "${Reads[0]}" | sed 's@.fastq.gz@@g')_fastqc.html

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 2 ]; do
		sleep 60s
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
#### Cutadapt  <a name="8"></a>
Primers were removed from the reads where present using Cutadapt.

NOTE:The reads are a mix of paired and single end samples.
```bash
for ReadDir in $(ls -d /data/users/theaven/nematode_project/raw_data/28S/*); do
	Task=CutAdapt
	ID=$(echo "$ReadDir" | cut -d '/' -f7,8 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
	OutDir="$(dirname "$ReadDir" | sed 's@raw_data@qc_data@g')/"$Task""
	Forward_Primer=ACAAGTACCGTGAGGGAAAGTTG
	Reverse_Primer=TCGGAAGGAACCAGCTACTA
	ExpectedOutput="$OutDir"/$(basename "${Reads[0]}" | sed 's@.fastq.gz@.trim.fastq.gz@g')

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 2 ]; do
		sleep 30s
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
	OutDir="$(dirname "$ReadDir" | sed 's@raw_data@qc_data@g')/"$Task""
	Forward_Primer=GCAAGTCTGGTGCCAGCAGC
	Reverse_Primer=CCGTGTTGAGTCAAATTAAG
	ExpectedOutput="$OutDir"/$(basename "${Reads[0]}" | sed 's@.fastq.gz@.trim.fastq.gz@g')

	Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
	while [ "$Jobs" -gt 2 ]; do
		sleep 30s
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
```
## Placeholder <a name="4"></a>

## Placeholder <a name="5"></a>

## Placeholder <a name="6"></a>
