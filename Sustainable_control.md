# Innovative Pest Sustainable control to reduce the impact of bark beetle and weevil attacks on alpine forests

All commands executed from /data/users/theaven/Ips_jam_project unless othewise specified.

Look for amplicons with the ITS primers in the Ips genome, https://primerdigital.com/tools/epcr.html returned no amplicons for the ITS primers.
```bash
echo ITS1-ITS4 TCCGTAGGTGAACCTGCGG TCCTCCGCTTATTGATATGC >> primers.txt
echo ITS1_Fus-ITS4 TCCGTTGGTGAACCAGCGG TCCTCCGCTTATTGATATGC >> primers.txt
echo ITS1_Mal-ITS4 TCTGTAGGTGAACCTGCAG TCCTCCGCTTATTGATATGC >> primers.txt
echo ITS1-ITS4_Pyt TCCGTAGGTGAACCTGCGG TCCTCCGCTTATTAATATGC >> primers.txt
echo ITS1_Fus-ITS4_Pyt TCCGTTGGTGAACCAGCGG TCCTCCGCTTATTAATATGC >> primers.txt
echo ITS1_Mal-ITS4_Pyt TCTGTAGGTGAACCTGCAG TCCTCCGCTTATTAATATGC >> primers.txt

module load anaconda3
conda activate emboss
primersearch \
  -seqall ~/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna \
  -infile primers.txt \
  -mismatchpercent 10 \
  -outfile primer_results.txt

module load samtools/1.19.2-gcc-13.3.0-a2yhwkt
awk -v GENOME="/home/clusterusers/theaven/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna" -f /home/clusterusers/theaven/git_repos/Scripts/unibz/primersearch_to_faidx.awk /home/clusterusers/theaven/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna.fai primer_results.txt > extract_cmds.sh
bash extract_cmds.sh > amplicons.fasta

#Fungal ITS sequences were downloaded from UNITE (DOI: 10.15156/BIO/3301229) [Accessed 1800202026], to this .fasta I will add the Ips sequences which have been predicted to be amplified by our ITS1/4 primmers
echo JADDUH010000185.1 >> temp_id.txt
echo JADDUH010000200.1 >> temp_id.txt
echo JADDUH010000208.1 >> temp_id.txt
echo JADDUH010000234.1 >> temp_id.txt
echo JADDUH010000076.1 >> temp_id.txt

module load apptainer/1.4.1-gcc-13.3.0-3coysxn
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py \
        --id_file temp_id.txt \
        --input /home/clusterusers/theaven/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna \
        --output temp.fa

cat temp.fa >> sh_general_release_dynamic_19.02.2025.fasta

#Create bedfile to define regions for depletion mode, these cover all of the hits in the Ips sequences + 250bp buffer (probably unnecessary for our use case)
echo JADDUH010000185.1 32622 35166 >> ITS_exclude.bed
echo JADDUH010000200.1 12958 15502 >> ITS_exclude.bed
echo JADDUH010000208.1 29499 32043 >> ITS_exclude.bed
echo JADDUH010000234.1 14676 17220 >> ITS_exclude.bed
echo JADDUH010000076.1 332025  399906 >> ITS_exclude.bed

cat sh_general_release_dynamic_19.02.2025.fasta | cut -f1 -d '|' >> ITS_all.fasta

awk '/^>/ {printf "%s%d\n", $0, ++i; next} {print}' ITS_all.fasta >> unite_ITS_all+.fasta 

grep '>' unite_ITS_all+.fasta | tail

>JADDUH010000185.1102138
>JADDUH010000200.1102139
>JADDUH010000208.1102140
>JADDUH010000234.1102141
>JADDUH010000076.1102142

srun -p bioagri  -c 8 --mem 32G --pty bash
module load anaconda3
conda activate vsearch

vsearch --derep_fulllength unite_its_19.02.2025.fasta \
  --output unite_its_19.02.2025_derep.fasta \
  --sizeout \
  --uc unite_derep.uc

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.95 \
  --centroids unite_its_19.02.2025_derep_95.fasta \
  --uc unite_95.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.90 \
  --centroids unite_its_19.02.2025_derep_90.fasta \
  --uc unite_90.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.85 \
  --centroids unite_its_19.02.2025_derep_85.fasta \
  --uc unite_85.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.80 \
  --centroids unite_its_19.02.2025_derep_80.fasta \
  --uc unite_80.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.75 \
  --centroids unite_its_19.02.2025_derep_75.fasta \
  --uc unite_75.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.70 \
  --centroids unite_its_19.02.2025_derep_70.fasta \
  --uc unite_70.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.65 \
  --centroids unite_its_19.02.2025_derep_65.fasta \
  --uc unite_65.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.60 \
  --centroids unite_its_19.02.2025_derep_60.fasta \
  --uc unite_60.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.55 \
  --centroids unite_its_19.02.2025_derep_55.fasta \
  --uc unite_55.uc \
  --threads 8

vsearch --cluster_fast unite_its_19.02.2025_derep.fasta \
  --id 0.50 \
  --centroids unite_its_19.02.2025_derep_50.fasta \
  --uc unite_50.uc \
  --threads 8

 for file in unite_its_19.02.2025_derep_*.fasta; do
  [ -e "$file" ] || continue

  echo "$file"
  grep -c '^>' "$file"

  out="${file%.fasta}_uni.fasta"

  awk 'BEGIN{i=0}
       /^>/{
         i++
         split($0, a, "|")     # keep only part before first |
         printf "%s_%d\n", a[1], i
         next
       }
       {print}' "$file" > "$out"
  cat temp.fa >> "$out"
done

unite_its_19.02.2025_derep_50.fasta
181

unite_its_19.02.2025_derep_55.fasta
647

unite_its_19.02.2025_derep_60.fasta
1662

unite_its_19.02.2025_derep_65.fasta
3576

unite_its_19.02.2025_derep_70.fasta
6241

unite_its_19.02.2025_derep_75.fasta
10023

unite_its_19.02.2025_derep_80.fasta
15439

unite_its_19.02.2025_derep_85.fasta
23587

unite_its_19.02.2025_derep_90.fasta
36270

unite_its_19.02.2025_derep_95.fasta
60303


```
```bash
#Compress and upload 16S files, check before and after upload
tar -cvzf "\\wsl.localhost\Ubuntu\home\tcheaven\16S.tar.gz" -C "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS" 16S
find /home/tcheaven/16S.tar.gz -type f -exec md5sum {} \; > /home/tcheaven/linux_16S_checksums.txt #9dd5f4045dcbbf4331a5103b2ecea0d1
find /data/users/theaven/Ips_jam_project/raw_data/minion/16S/16S.tar.gz -type f -exec md5sum {} \; > /data/users/theaven/Ips_jam_project/raw_data/minion/16S/linux_16S_checksums.txt #9dd5f4045dcbbf4331a5103b2ecea0d1
#checksums of the .tar.gz files are the same before and after upload
tar -xzvf 16S.tar.gz

##################################################################################################

#upload of the compressed ITS file keeps failing, presumably because it is too large, uploading individula directories instead.
#Check key files:
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24\20260225-TOMH-IpsITS-1-24\20260225_1922_MN41812_FBF74689_a36ad5b2\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-2\20260225-TOMH-IpsITS-1-24-2\20260225_1944_MN41812_FBF74689_913d586f\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums2.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-3\20260225-TOMH-IpsITS-1-24-3\20260225_2000_MN41812_FBF74689_a1aa9534\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums3.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-4\20260225-TOMH-IpsITS-1-24-4\20260225_2015_MN41812_FBF74689_8467c4fc\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums4.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-7\20260225-TOMH-IpsITS-1-24-7\20260225_2052_MN41812_FBF74689_ae8e4602\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums7.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-8\20260225-TOMH-IpsITS-1-24-8\20260225_2115_MN41812_FBF74689_51d6a6cc\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums8.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-9\20260225-TOMH-IpsITS-1-24-9\20260225_2129_MN41812_FBF74689_7717758a\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums9.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-10\20260225-TOMH-IpsITS-1-24-10\20260226_1132_MN41812_FBF74689_37324c24\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums10.txt"
Get-FileHash "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS\ITS\20260225-TOMH-IpsITS-1-24-12\20260225-TOMH-IpsITS-1-24-12\20260226_2313_MN41812_FBF74689_c024ebc9\pod5\*" -Algorithm MD5 |
ForEach-Object { "$($_.Hash)  $(Split-Path $_.Path -Leaf)" } |
Out-File "C:\Users\THeaven\windows_ITS_checksums12.txt"

find /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/*/*/*/pod5/*.pod5 -type f -exec md5sum {} \; >> /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/linux_ITS_checksums.txt
#checksums of the .pod5 files are the same before and after upload
```
#### Basecalling

```bash
mkdir /data/users/theaven/Ips_jam_project/raw_data/minion/16S/pod5
mkdir /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/pod5

ln -s /data/users/theaven/Ips_jam_project/raw_data/minion/16S/*/*/*/pod5/*.pod5 /data/users/theaven/Ips_jam_project/raw_data/minion/16S/pod5/.
ln -s /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/*/*/*/pod5/*.pod5 /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/pod5/.

tar -cvzf "\\wsl.localhost\Ubuntu\home\tcheaven\ITS2.tar.gz" -C "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\nanopore_16S_ITS" ITS


srun  -p gpu-low  -c 16 --gres=gpu:1 --mem=64G --account=shame --pty bash
module load apptainer/1.4.1-gcc-13.3.0-3                         
module load nextflow/23.10.1-gcc-12.1.0
nextflow run epi2me-labs/wf-basecalling --help
nextflow pull epi2me-labs/wf-basecalling
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-basecalling/wf-basecalling-demo.tar.gz
tar -xzvf wf-basecalling-demo.tar.gz

mkdir -p /data/users/theaven/Ips_jam_project/.singularity-cache
mkdir -p /data/users/theaven/Ips_jam_project/.apptainer-cache
export NXF_SINGULARITY_CACHEDIR=/data/users/theaven/Ips_jam_project/.singularity-cache
export APPTAINER_CACHEDIR=/data/users/theaven/Ips_jam_project/.apptainer-cache
nextflow run epi2me-labs/wf-basecalling \
  --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
  --dorado_ext 'pod5' \
  --input 'wf-basecalling-demo/input' \
  --ref 'wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
  --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
  -profile singularity \
  -resume

nextflow run epi2me-labs/wf-basecalling \
  --max_cpus 16 --ubam_map_threads 16 --ubam_sort_threads 8 --ubam_bam2fq_threads 4 --merge_threads 4 \
  --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0' \
  --dorado_ext 'pod5' \
  --cuda_device 'cuda:all' \
  --input '/data/users/theaven/Ips_jam_project/raw_data/minion/16S/pod5' \
  --barcode_kit 'SQK-MAB114-24' \
  -profile singularity \
  --out_dir /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls --output_fmt fastq \
  -resume

screen -S dorado
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/16S/pod5 /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/pod5); do
  Task=Dorado
  InDir="$Dir"
  OutDir=$(dirname $Dir)/basecalls
  OutFmt=fastq
  Barcode=SQK-MAB114-24
  Modification_model=NA
  ExpectedOutput="$OutDir"/out.fastq

  Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  while [ "$Jobs" -gt 0 ]; do
    sleep 600s
    printf "."
    Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  done

  if [ ! -s "$ExpectedOutput" ]; then
    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_dorado.sh "$InDir" "$OutDir" "$OutFmt" "$Barcode" "$Modification_model")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

#Combine the different runs into single files per demultiplexed sample:
for dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/demuxed/*/*/*/fastq_pass/*); do
  barcode=$(basename $dir)
  cat "$dir"/*.fastq >> /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/"$barcode".fastq
done

for dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/demuxed/*/*/*/fastq_pass/*); do
  barcode=$(basename $dir)
  cat "$dir"/*.fastq >> /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/"$barcode".fastq
done

```
#### EPI2ME - wf-16s

https://epi2me.nanoporetech.com/epi2me-docs/workflows/wf-16s/#analysing-its-amplicons

```bash
#Test the pipeline
screen -S wf-16s
srun  -p bioagri  -c 12 --mem=32G --account=shame --pty bash
module load apptainer/1.4.1-gcc-13.3.0-3
module load nextflow/23.10.1-gcc-12.1.0


nextflow run epi2me-labs/wf-16s --help
nextflow pull epi2me-labs/wf-16s
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-16s/wf-16s-demo.tar.gz
tar -xzvf wf-16s-demo.tar.gz

nextflow run epi2me-labs/wf-16s \
    --fastq 'wf-16s-demo/test_data' \
    --minimap2_by_reference \
    -profile singularity

#Make expected file structure for multiplexed reads
for dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/*/basecalls/demuxed/*/*/*/fastq_pass/*); do
  barcode=$(basename $dir)
  Out=$(echo "$dir" | cut -d '/' -f1,2,3,4,5,6,7,8 | sed 's@raw_data/minion@wf-16s/in@g')/"$barcode"
  mkdir -p "$Out"
  ln -s "$dir"/* "$Out"/.
done

echo barcode,alias > 20260323_sample_sheet.csv
echo barcode01,SHBB01881-1 >> 20260323_sample_sheet.csv
echo barcode02,SHBB01881-2 >> 20260323_sample_sheet.csv
echo barcode03,SHBB01881-3 >> 20260323_sample_sheet.csv
echo barcode04,SHBB01882-1 >> 20260323_sample_sheet.csv
echo barcode05,SHBB01882-2 >> 20260323_sample_sheet.csv
echo barcode06,SHBB01882-4 >> 20260323_sample_sheet.csv
echo barcode07,SHBB01884-3 >> 20260323_sample_sheet.csv
echo barcode08,SHBB01887-2 >> 20260323_sample_sheet.csv
echo barcode09,SHBB01888-1 >> 20260323_sample_sheet.csv
echo barcode10,SHBB01888-3 >> 20260323_sample_sheet.csv
echo barcode11,SHBB01888-4 >> 20260323_sample_sheet.csv
echo barcode12,SHBB01889-1 >> 20260323_sample_sheet.csv
echo barcode13,SHBB01890-1 >> 20260323_sample_sheet.csv
echo barcode14,SHBB01891-1 >> 20260323_sample_sheet.csv
echo barcode15,SHBB01895-1 >> 20260323_sample_sheet.csv
echo barcode16,SHBB01898-1 >> 20260323_sample_sheet.csv
echo barcode17,SHBB01899-1 >> 20260323_sample_sheet.csv
echo barcode18,SHBB01900-1 >> 20260323_sample_sheet.csv
echo barcode19,SHBB01903-1 >> 20260323_sample_sheet.csv
echo barcode20,SHBB01909-1 >> 20260323_sample_sheet.csv
echo barcode21,SHBB01910-1 >> 20260323_sample_sheet.csv
echo barcode22,SHBB01911-1 >> 20260323_sample_sheet.csv
echo barcode23,SHBB01914-1 >> 20260323_sample_sheet.csv
echo barcode24,SHBB01915-1 >> 20260323_sample_sheet.csv

#Run the pipeline

#16S with all reads, minimap
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/16S); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/16S/inclusive
  Classifier=minimap2
  Database=ncbi_16s_18s
  Exclude=N
  Max_len=100000
  Min_len=10
  Abundance=0.0
  Unclassified=true
  ExpectedOutput="$OutDir"/out.fastq

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EPI2ME_wf-16s.sh "$InDir" "$Sample_Sheet" "$Amplicon" "$OutDir" "$Classifier" --database "$Database" --exclude "$Exclude" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance" --unclassified "$Unclassified")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

#ITS with all reads, minimap
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/ITS); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/ITS/inclusive
  Classifier=minimap2
  Database=ncbi_16s_18s_28s_ITS
  Exclude=N
  Max_len=100000
  Min_len=10
  Abundance=0.0
  Unclassified=true
  ExpectedOutput="$OutDir"/out.fastq

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EPI2ME_wf-16s.sh "$InDir" "$Sample_Sheet" "$Amplicon" "$OutDir" "$Classifier" --database "$Database" --exclude "$Exclude" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance" --unclassified "$Unclassified")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

#16S with expected size reads, removing host, minimap
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/16S); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/16S/exclusive
  Classifier=minimap2
  Database=ncbi_16s_18s
  Exclude=/home/clusterusers/theaven/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna
  Max_len=2000
  Min_len=1000
  Abundance=0.01
  Unclassified=false
  ExpectedOutput="$OutDir"/out.fastq

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EPI2ME_wf-16s.sh "$InDir" "$Sample_Sheet" "$Amplicon" "$OutDir" "$Classifier" --database "$Database" --exclude "$Exclude" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance" --unclassified "$Unclassified")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

#ITS with expected size reads, removing host, minimap
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/ITS); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/ITS/exclusive
  Classifier=minimap2
  Database=ncbi_16s_18s_28s_ITS
  Exclude=/home/clusterusers/theaven/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna
  Max_len=1000
  Min_len=300
  Abundance=0.01
  Unclassified=false
  ExpectedOutput="$OutDir"/out.fastq

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EPI2ME_wf-16s.sh "$InDir" "$Sample_Sheet" "$Amplicon" "$OutDir" "$Classifier" --database "$Database" --exclude "$Exclude" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance" --unclassified "$Unclassified")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

#16S with expected size reads, removing host, kraken
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/16S); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/16S/exclusive-kraken
  Classifier=kraken2
  Database=ncbi_16s_18s
  Exclude=/home/clusterusers/theaven/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna
  Max_len=2000
  Min_len=1000
  Abundance=0.01
  Unclassified=false
  ExpectedOutput="$OutDir"/out.fastq

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EPI2ME_wf-16s.sh "$InDir" "$Sample_Sheet" "$Amplicon" "$OutDir" "$Classifier" --database "$Database" --exclude "$Exclude" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance" --unclassified "$Unclassified")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

#ITS with expected size reads, removing host, kraken
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/ITS); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/ITS/exclusive-kraken
  Classifier=kraken2
  Database=ncbi_16s_18s_28s_ITS
  Exclude=/home/clusterusers/theaven/genomes/Ips/typographus/GCA_016097725.1/GCA_016097725.1_CZU_Ityp_1.0_genomic.fna
  Max_len=1000
  Min_len=300
  Abundance=0.01
  Unclassified=false
  ExpectedOutput="$OutDir"/out.fastq

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EPI2ME_wf-16s.sh "$InDir" "$Sample_Sheet" "$Amplicon" "$OutDir" "$Classifier" --database "$Database" --exclude "$Exclude" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance" --unclassified "$Unclassified")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done
```