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

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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

Sequencing was run with fast basecalling, raw .POD5 files were output which we will now use for basecalling with the highest accuracy settings with dorado. Dorado is designed by Oxford Nanopore and replaces guppy, it performs basecalling, demultiplexing, and also adapter trimming by default - there is therefore no need to use secondary adapter trimmers, eg. porechop.

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
#### Primer trimming

```bash
module load anaconda3
conda activate seqkit-2.10
seqkit stats /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/*.fastq
```
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode01.fastq     FASTQ   DNA    111,241  149,580,502        8  1,344.7   12,844
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode02.fastq     FASTQ   DNA    140,325  179,014,046       11  1,275.7   10,432
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode03.fastq     FASTQ   DNA    110,977  148,984,323        8  1,342.5   17,140
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode04.fastq     FASTQ   DNA    113,067  152,202,613       16  1,346.1    7,277
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode05.fastq     FASTQ   DNA     30,957   41,353,842       27  1,335.8    5,144
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode06.fastq     FASTQ   DNA     41,886   57,011,854       16  1,361.1   10,744
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode07.fastq     FASTQ   DNA      5,965    8,044,267       12  1,348.6    3,105
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode08.fastq     FASTQ   DNA    102,661  138,795,354       19    1,352   10,865
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode09.fastq     FASTQ   DNA     79,755  104,961,509       16    1,316    4,483
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode10.fastq     FASTQ   DNA    105,901  142,765,821        8  1,348.1    4,583
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode11.fastq     FASTQ   DNA     20,363   27,084,171       15  1,330.1   13,363
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode12.fastq     FASTQ   DNA    112,612  151,234,068       19    1,343   12,956
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode13.fastq     FASTQ   DNA     98,895  130,059,284       18  1,315.1    8,139
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode14.fastq     FASTQ   DNA    107,198  145,827,884       20  1,360.4   14,618
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode15.fastq     FASTQ   DNA    105,822  144,867,178       15    1,369   18,744
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode16.fastq     FASTQ   DNA     30,762   40,438,892       14  1,314.6   14,102
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode17.fastq     FASTQ   DNA    199,358  261,826,410       16  1,313.3   24,287
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode18.fastq     FASTQ   DNA    128,432  170,680,199       10    1,329    6,077
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode19.fastq     FASTQ   DNA    241,115  325,755,979       11    1,351   15,792
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode20.fastq     FASTQ   DNA    160,949  220,457,016       38  1,369.7   14,937
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode21.fastq     FASTQ   DNA    310,507  401,804,911        7    1,294   15,106
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode22.fastq     FASTQ   DNA    229,969  304,474,323       13    1,324   11,370
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode23.fastq     FASTQ   DNA    183,902  252,614,214        7  1,373.6   13,551
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode24.fastq     FASTQ   DNA     67,645   97,286,004       16  1,438.2   19,715
/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/unclassified.fastq  FASTQ   DNA    213,830  319,071,976        5  1,492.2   27,940

```bash
for Reads in /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/*.fastq; do
  Task=CutAdapt
  OutDir="$(dirname "$Reads" | sed 's@raw_data@qc_data@g')"/"$Task"
  Amplicon=16S
  ExpectedOutput="$OutDir"/$(basename "$Reads" | sed 's@.fastq@.trim.fastq@g')

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p "$OutDir"
    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_cutadapt_nanopore_16S-ITS.sh "$OutDir" "$Amplicon" "$Reads")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for file in $(ls /home/clusterusers/theaven/slurm_records/slurm.205*.out); do
grep 'ReadFiles: ' $file >> /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/cutadapt_report.txt
grep -A 12 '=== Summary ===' $file >> /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/cutadapt_report.txt
echo -e '\n\n' >> /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/cutadapt_report.txt
done

#62.4% of reads dropped, 8.18% because of no primer sequence in the reads, 0.29% too long, 51.73% tool long on average

seqkit stats /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/*.fastq
```
file                                                                                               format  type  num_seqs      sum_len  min_len  avg_len  max_len
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode01.trim.fastq     FASTQ   DNA     16,978   24,076,817      800  1,418.1    1,971
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode02.trim.fastq     FASTQ   DNA     27,600   38,676,559      800  1,401.3    1,988
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode03.trim.fastq     FASTQ   DNA     19,959   28,124,779      800  1,409.1    1,937
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode04.trim.fastq     FASTQ   DNA     24,447   34,230,304      800  1,400.2    1,989
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode05.trim.fastq     FASTQ   DNA     15,997   22,425,799      800  1,401.9    1,985
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode06.trim.fastq     FASTQ   DNA     21,036   29,927,711      800  1,422.7    1,946
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode07.trim.fastq     FASTQ   DNA      3,660    5,234,820      801  1,430.3    1,824
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode08.trim.fastq     FASTQ   DNA     49,220   68,972,949      800  1,401.3    1,988
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode09.trim.fastq     FASTQ   DNA     32,593   46,190,733      800  1,417.2    1,992
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode10.trim.fastq     FASTQ   DNA     45,121   64,443,143      800  1,428.2    1,936
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode11.trim.fastq     FASTQ   DNA      9,511   13,334,504      800    1,402    1,999
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode12.trim.fastq     FASTQ   DNA     49,949   70,831,260      800  1,418.1    1,909
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode13.trim.fastq     FASTQ   DNA     45,020   63,098,355      800  1,401.6    1,989
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode14.trim.fastq     FASTQ   DNA     48,393   69,326,418      800  1,432.6    1,992
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode15.trim.fastq     FASTQ   DNA     47,294   67,801,961      800  1,433.6    1,979
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode16.trim.fastq     FASTQ   DNA     13,582   18,954,530      801  1,395.6    1,966
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode17.trim.fastq     FASTQ   DNA     77,561  108,733,977      800  1,401.9    1,986
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode18.trim.fastq     FASTQ   DNA     53,409   75,564,777      800  1,414.8    1,994
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode19.trim.fastq     FASTQ   DNA    100,699  143,511,771      800  1,425.2    1,989
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode20.trim.fastq     FASTQ   DNA     63,954   90,115,020      800  1,409.1    1,988
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode21.trim.fastq     FASTQ   DNA    127,780  176,686,645      800  1,382.7    1,998
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode22.trim.fastq     FASTQ   DNA     87,432  122,411,101      800  1,400.1    1,991
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode23.trim.fastq     FASTQ   DNA     68,561   95,610,238      800  1,394.5    1,982
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode24.trim.fastq     FASTQ   DNA     24,936   35,577,616      800  1,426.8    1,979
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/unclassified.trim.fastq  FASTQ   DNA     74,297  106,173,199      800    1,429    2,000

```bash
seqkit stats /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/*.fastq
```
file                                                                                  format  type  num_seqs      sum_len  min_len  avg_len  max_len
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode01.fastq     FASTQ   DNA     80,575   73,348,711       20    910.3   11,814
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode02.fastq     FASTQ   DNA      7,560    5,247,563       16    694.1    4,979
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode03.fastq     FASTQ   DNA    214,544  131,426,422       20    612.6    8,757
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode04.fastq     FASTQ   DNA     73,598   56,328,797       16    765.4   11,084
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode05.fastq     FASTQ   DNA     55,288   39,793,018       20    719.7    6,221
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode06.fastq     FASTQ   DNA    121,694   78,177,092       14    642.4    8,733
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode07.fastq     FASTQ   DNA    101,232   71,750,291       21    708.8   13,414
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode08.fastq     FASTQ   DNA     79,365   82,466,194       18  1,039.1   15,927
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode09.fastq     FASTQ   DNA     10,475    8,866,068        8    846.4   10,792
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode10.fastq     FASTQ   DNA    105,718  114,033,817       16  1,078.7    4,142
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode11.fastq     FASTQ   DNA    132,047  115,624,945       15    875.6    4,165
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode12.fastq     FASTQ   DNA     37,096   50,004,586       14    1,348   11,690
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode13.fastq     FASTQ   DNA     31,651   23,918,465       21    755.7    7,771
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode14.fastq     FASTQ   DNA     73,219  104,730,913       28  1,430.4   15,737
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode15.fastq     FASTQ   DNA    143,152  157,399,347       20  1,099.5   11,807
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode16.fastq     FASTQ   DNA     23,665   32,631,787       20  1,378.9   21,641
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode17.fastq     FASTQ   DNA     88,673   79,828,607       31    900.3    5,803
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode18.fastq     FASTQ   DNA     58,003   52,526,223       22    905.6    5,624
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode19.fastq     FASTQ   DNA    122,800  110,476,251       31    899.6   10,503
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode20.fastq     FASTQ   DNA     74,967   69,095,985       20    921.7    7,620
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode21.fastq     FASTQ   DNA    115,329  105,720,604       19    916.7   18,451
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode22.fastq     FASTQ   DNA     35,859   37,231,366       12  1,038.3    5,341
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode23.fastq     FASTQ   DNA    207,039  151,490,589       23    731.7   10,753
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/barcode24.fastq     FASTQ   DNA    179,113  171,499,626       21    957.5   10,711
/data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/unclassified.fastq  FASTQ   DNA     88,567   96,722,643       65  1,092.1  123,248

```bash
for Reads in /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/*.fastq; do
  Task=CutAdapt
  OutDir="$(dirname "$Reads" | sed 's@raw_data@qc_data@g')"/"$Task"
  Amplicon=ITS
  ExpectedOutput="$OutDir"/$(basename "$Reads" | sed 's@.fastq@.trim.fastq@g')

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p "$OutDir"
    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_cutadapt_nanopore_16S-ITS.sh "$OutDir" "$Amplicon" "$Reads")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for file in $(ls /home/clusterusers/theaven/slurm_records/slurm.205*.out); do
grep 'ReadFiles: ' $file >> /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/cutadapt_report.txt
grep -A 12 '=== Summary ===' $file >> /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/cutadapt_report.txt
echo -e '\n\n' >> /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/cutadapt_report.txt
done

#70.4% of reads dropped, 1.65% because of no primer sequence in the reads, 7.26% too long, 59.9% tool long on average

seqkit stats /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/*.fastq
```
file                                                                                               format  type  num_seqs     sum_len  min_len  avg_len  max_len
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode01.trim.fastq     FASTQ   DNA     20,965  12,666,280      300    604.2    1,495
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode02.trim.fastq     FASTQ   DNA      5,253   3,033,425      300    577.5    1,426
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode03.trim.fastq     FASTQ   DNA     74,826  42,646,806      300    569.9    1,491
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode04.trim.fastq     FASTQ   DNA     25,306  15,886,093      300    627.8    1,499
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode05.trim.fastq     FASTQ   DNA     22,846  13,930,136      300    609.7    1,499
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode06.trim.fastq     FASTQ   DNA     52,200  31,370,661      300      601    1,491
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode07.trim.fastq     FASTQ   DNA     30,256  17,299,892      300    571.8    1,485
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode08.trim.fastq     FASTQ   DNA     19,459  12,875,250      300    661.7    1,497
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode09.trim.fastq     FASTQ   DNA      5,318   4,010,328      300    754.1    1,428
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode10.trim.fastq     FASTQ   DNA     30,588  24,411,630      300    798.1    1,495
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode11.trim.fastq     FASTQ   DNA     45,968  29,738,074      300    646.9    1,498
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode12.trim.fastq     FASTQ   DNA      5,874   3,384,437      300    576.2    1,494
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode13.trim.fastq     FASTQ   DNA     12,595   7,968,252      300    632.7    1,478
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode14.trim.fastq     FASTQ   DNA     10,377   5,743,763      300    553.5    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode15.trim.fastq     FASTQ   DNA     32,177  18,797,165      300    584.2    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode16.trim.fastq     FASTQ   DNA      4,604   3,790,759      300    823.4    1,438
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode17.trim.fastq     FASTQ   DNA     28,388  17,546,669      300    618.1    1,499
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode18.trim.fastq     FASTQ   DNA     12,542  10,463,073      300    834.2    1,492
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode19.trim.fastq     FASTQ   DNA     29,290  16,506,781      300    563.6    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode20.trim.fastq     FASTQ   DNA     24,153  18,410,205      300    762.2    1,493
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode21.trim.fastq     FASTQ   DNA     37,020  26,769,559      300    723.1    1,498
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode22.trim.fastq     FASTQ   DNA      8,182   4,857,078      300    593.6    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode23.trim.fastq     FASTQ   DNA     63,012  40,311,448      300    639.7    1,478
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/barcode24.trim.fastq     FASTQ   DNA     44,562  26,405,739      300    592.6    1,497
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/unclassified.trim.fastq  FASTQ   DNA     23,231  15,720,375      300    676.7    1,500

#### Filtlong

```bash
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/*.fastq; do
  Task=Filtlong
  OutDir="$(dirname "$Reads")"/"$Task"
  Amplicon=16S
  ExpectedOutput="$OutDir"/$(basename "$Reads" | sed 's@.fastq@.trim.fastq@g')

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p "$OutDir"
    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_filtlong_nanopore_16S-ITS.sh "$OutDir" "$Amplicon" "$Reads")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

seqkit stats /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/*.fastq.gz
```
file                                                                                                                    format  type  num_seqs      sum_len  min_len  avg_len  max_len
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode01.trim.filtlong.fastq.gz     FASTQ   DNA     16,978   24,076,817      800  1,418.1    1,971
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode02.trim.filtlong.fastq.gz     FASTQ   DNA     27,600   38,676,559      800  1,401.3    1,988
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode03.trim.filtlong.fastq.gz     FASTQ   DNA     19,959   28,124,779      800  1,409.1    1,937
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode04.trim.filtlong.fastq.gz     FASTQ   DNA     24,447   34,230,304      800  1,400.2    1,989
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode05.trim.filtlong.fastq.gz     FASTQ   DNA     15,997   22,425,799      800  1,401.9    1,985
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode06.trim.filtlong.fastq.gz     FASTQ   DNA     21,036   29,927,711      800  1,422.7    1,946
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode07.trim.filtlong.fastq.gz     FASTQ   DNA      3,660    5,234,820      801  1,430.3    1,824
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode08.trim.filtlong.fastq.gz     FASTQ   DNA     49,220   68,972,949      800  1,401.3    1,988
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode09.trim.filtlong.fastq.gz     FASTQ   DNA     32,593   46,190,733      800  1,417.2    1,992
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode10.trim.filtlong.fastq.gz     FASTQ   DNA     45,121   64,443,143      800  1,428.2    1,936
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode11.trim.filtlong.fastq.gz     FASTQ   DNA      9,511   13,334,504      800    1,402    1,999
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode12.trim.filtlong.fastq.gz     FASTQ   DNA     49,949   70,831,260      800  1,418.1    1,909
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode13.trim.filtlong.fastq.gz     FASTQ   DNA     45,020   63,098,355      800  1,401.6    1,989
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode14.trim.filtlong.fastq.gz     FASTQ   DNA     48,393   69,326,418      800  1,432.6    1,992
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode15.trim.filtlong.fastq.gz     FASTQ   DNA     47,294   67,801,961      800  1,433.6    1,979
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode16.trim.filtlong.fastq.gz     FASTQ   DNA     13,582   18,954,530      801  1,395.6    1,966
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode17.trim.filtlong.fastq.gz     FASTQ   DNA     77,561  108,733,977      800  1,401.9    1,986
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode18.trim.filtlong.fastq.gz     FASTQ   DNA     53,409   75,564,777      800  1,414.8    1,994
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode19.trim.filtlong.fastq.gz     FASTQ   DNA    100,699  143,511,771      800  1,425.2    1,989
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode20.trim.filtlong.fastq.gz     FASTQ   DNA     63,954   90,115,020      800  1,409.1    1,988
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode21.trim.filtlong.fastq.gz     FASTQ   DNA    127,780  176,686,645      800  1,382.7    1,998
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode22.trim.filtlong.fastq.gz     FASTQ   DNA     87,432  122,411,101      800  1,400.1    1,991
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode23.trim.filtlong.fastq.gz     FASTQ   DNA     68,561   95,610,238      800  1,394.5    1,982
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/barcode24.trim.filtlong.fastq.gz     FASTQ   DNA     24,936   35,577,616      800  1,426.8    1,979
/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/unclassified.trim.filtlong.fastq.gz  FASTQ   DNA     74,297  106,173,199      800    1,429    2,000

```bash
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/*.fastq; do
  Task=Filtlong
  OutDir="$(dirname "$Reads")"/"$Task"
  Amplicon=ITS
  ExpectedOutput="$OutDir"/$(basename "$Reads" | sed 's@.fastq@.trim.fastq@g')

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p "$OutDir"
    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_filtlong_nanopore_16S-ITS.sh "$OutDir" "$Amplicon" "$Reads")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

seqkit stats /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.fastq.gz
```
file                                                                                                                    format  type  num_seqs     sum_len  min_len  avg_len  max_len
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode01.trim.filtlong.fastq.gz     FASTQ   DNA     20,965  12,666,280      300    604.2    1,495
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode02.trim.filtlong.fastq.gz     FASTQ   DNA      5,253   3,033,425      300    577.5    1,426
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode03.trim.filtlong.fastq.gz     FASTQ   DNA     74,826  42,646,806      300    569.9    1,491
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode04.trim.filtlong.fastq.gz     FASTQ   DNA     25,306  15,886,093      300    627.8    1,499
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode05.trim.filtlong.fastq.gz     FASTQ   DNA     22,846  13,930,136      300    609.7    1,499
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode06.trim.filtlong.fastq.gz     FASTQ   DNA     52,200  31,370,661      300      601    1,491
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode07.trim.filtlong.fastq.gz     FASTQ   DNA     30,256  17,299,892      300    571.8    1,485
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode08.trim.filtlong.fastq.gz     FASTQ   DNA     19,459  12,875,250      300    661.7    1,497
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode09.trim.filtlong.fastq.gz     FASTQ   DNA      5,318   4,010,328      300    754.1    1,428
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode10.trim.filtlong.fastq.gz     FASTQ   DNA     30,588  24,411,630      300    798.1    1,495
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode11.trim.filtlong.fastq.gz     FASTQ   DNA     45,968  29,738,074      300    646.9    1,498
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode12.trim.filtlong.fastq.gz     FASTQ   DNA      5,874   3,384,437      300    576.2    1,494
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode13.trim.filtlong.fastq.gz     FASTQ   DNA     12,595   7,968,252      300    632.7    1,478
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode14.trim.filtlong.fastq.gz     FASTQ   DNA     10,377   5,743,763      300    553.5    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode15.trim.filtlong.fastq.gz     FASTQ   DNA     32,177  18,797,165      300    584.2    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode16.trim.filtlong.fastq.gz     FASTQ   DNA      4,604   3,790,759      300    823.4    1,438
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode17.trim.filtlong.fastq.gz     FASTQ   DNA     28,388  17,546,669      300    618.1    1,499
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode18.trim.filtlong.fastq.gz     FASTQ   DNA     12,542  10,463,073      300    834.2    1,492
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode19.trim.filtlong.fastq.gz     FASTQ   DNA     29,290  16,506,781      300    563.6    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode20.trim.filtlong.fastq.gz     FASTQ   DNA     24,153  18,410,205      300    762.2    1,493
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode21.trim.filtlong.fastq.gz     FASTQ   DNA     37,020  26,769,559      300    723.1    1,498
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode22.trim.filtlong.fastq.gz     FASTQ   DNA      8,182   4,857,078      300    593.6    1,500
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode23.trim.filtlong.fastq.gz     FASTQ   DNA     63,012  40,311,448      300    639.7    1,478
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/barcode24.trim.filtlong.fastq.gz     FASTQ   DNA     44,562  26,405,739      300    592.6    1,497
/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/unclassified.trim.filtlong.fastq.gz  FASTQ   DNA     23,060  15,657,253      300      679    1,500


The only reads dropped by the --min_mean_q 70 quality threshold of Filtlong are from the unclassified file.

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

Seperate with and without depletion
```bash
#Make expected file structure for multiplexed reads
for dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/*/basecalls/demuxed/*-9/*/*/fastq_pass/*); do
  barcode=$(basename $dir)
  Out=$(echo "$dir" | cut -d '/' -f1,2,3,4,5,6,7,8 | sed 's@raw_data/minion/ITS@wf-16s/in/ITS-nodep@g')/"$barcode"
  mkdir -p "$Out"
  ln -s "$dir"/* "$Out"/.
done

for dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/*/basecalls/demuxed/*-10/*/*/fastq_pass/*); do
  barcode=$(basename $dir)
  Out=$(echo "$dir" | cut -d '/' -f1,2,3,4,5,6,7,8 | sed 's@raw_data/minion/ITS@wf-16s/in/ITS-dep@g')/"$barcode"
  mkdir -p "$Out"
  ln -s "$dir"/* "$Out"/.
done

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/ITS-dep); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/ITS/exclusive-kraken-dep
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

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/ITS-nodep); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/ITS/exclusive-kraken-nodep
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
```bash
mkdir /data/users/theaven/download_20260330
for file in $(ls /data/users/theaven/Ips_jam_project/wf-16s/out/*/*/*-report.html); do
cp $file /data/users/theaven/download_20260330/$(echo $file | sed 's@/data/users/theaven/download_20260330@@g' | sed 's@/@-@g' )
done
```
***Repeat EPI2ME with filtered reads***
```bash
#Make expected file structure for multiplexed reads
for dir in $(ls /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/*.trim.fastq); do
  barcode=$(basename $dir | sed 's@.trim.fastq@@g')
  Out=$(echo "$dir" | cut -d '/' -f1,2,3,4,5,6,7,8 | sed 's@qc_data/minion/16S@wf-16s/in/16S-trim@g')/"$barcode"
  mkdir -p "$Out"
  ln -s "$dir" "$Out"
done

for dir in $(ls /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/*.trim.fastq); do
  barcode=$(basename $dir | sed 's@.trim.fastq@@g')
  Out=$(echo "$dir" | cut -d '/' -f1,2,3,4,5,6,7,8 | sed 's@qc_data/minion/ITS@wf-16s/in/ITS-trim@g')/"$barcode"
  mkdir -p "$Out"
  ln -s "$dir" "$Out"
done

#16S with expected size reads, removing host, minimap
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/16S-trim); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/16S-trim/exclusive
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
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/ITS-trim); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/ITS-trim/exclusive
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
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/16S-trim); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/16S-trim/exclusive-kraken
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
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/wf-16s/in/ITS-trim); do
  Task=wf-16s
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  OutDir=/data/users/theaven/Ips_jam_project/wf-16s/out/ITS-trim/exclusive-kraken
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

## EMU

Emu is a relative abundance estimator for 16S genomic sequences. The method is optimized for error-prone full-length reads, but can also be utilized for short-read data.

Download databases, from https://osf.io/32sh5/overview and https://osf.io/56uf7/overview :

```bash
#Updated emu rrnDB v5.10 and NCBI March 2026 courtesy of @UFDuttonLab [01/04/2026]


#EMU prebuilt database contains a combination of rrnDB v5.6 and NCBI 16S RefSeq from 17 September, 2020. Taxonomy is also from NCBI on the same date. The resulting database contains 49,301 sequences from 17,555 unique bacterial and archaeal species. [01/02/2023]
mkdir /data/users/theaven/db/emu/emu
export EMU_DATABASE_DIR=/data/users/theaven/db/emu/emu
cd "$EMU_DATABASE_DIR"
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
tar -xvf emu.tar

#Updated emu rrnDB v5.10 and NCBI March 2026 courtesy of @UFDuttonLab. [01/04/2026]
mkdir /data/users/theaven/db/emu/emu2026
cd ../emu2026
osf -p 32sh5 fetch osfstorage/overview/species_taxid.fasta
osf -p 32sh5 fetch osfstorage/overview/taxonomy.tsv

#SILVA v138.1 has been pre-built for Emu v3.0+ from the DADA2 SILVA species-level database. [01/02/2023]
mkdir /data/users/theaven/db/emu/silva-138.1
cd ../silva-138.1
osf -p 56uf7 fetch osfstorage/emu-prebuilt/silva.tar
tar -xvf silva.tar

#SILVA v138.2 has been built courtesy of @maddne. This database has not yet been tested or validated with Emu. [16/01/2026]
mkdir /data/users/theaven/db/emu/silva
cd ../silva
osf -p 56uf7 fetch osfstorage/emu-prebuilt/silva-138.2.tar
tar -xvf silva-138.2.tar

#UNITE general fasta v8.3 fungi has been pre-built for Emu v3.0+. This database has not yet been tested or validated with Emu. [01/02/2023]
mkdir /data/users/theaven/db/emu/unite-fungi
cd ../unite-fungi
osf -p 56uf7 fetch osfstorage/emu-prebuilt/unite-fungi.tar
tar -xvf unite-fungi.tar

#UNITE general fasta v8.3 all eukaryotes has been pre-built for Emu v3.0+. This database has not yet been tested or validated with Emu. [01/02/2023]
mkdir /data/users/theaven/db/emu/unite-all
cd ../unite-all
osf -p 56uf7 fetch osfstorage/emu-prebuilt/unite-all.tar
tar -xvf unite-all.tar
```
```bash
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  Type=map-ont
  Database=/data/users/theaven/db/emu/emu2026
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  Type=map-ont
  Database=/data/users/theaven/db/emu/emu
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  Type=map-ont
  Database=/data/users/theaven/db/emu/silva
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=16S
  Type=map-ont
  Database=/data/users/theaven/db/emu/silva-138.1
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

###################################################################################################################

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  Type=map-ont
  Database=/data/users/theaven/db/emu/silva
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  Type=map-ont
  Database=/data/users/theaven/db/emu/silva-138.1
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  Type=map-ont
  Database=/data/users/theaven/db/emu/unite-all
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet.csv
  Amplicon=ITS
  Type=map-ont
  Database=/data/users/theaven/db/emu/unite-fungi
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/1
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done
```

"unmapped" reads are reads that did not result in a mapping to the provided database with minimap2. "Filtered mapped" reads are those that were mapped with minimap2, but all alignments for the given query (read) were filtered via the align-len and percent identity (pid) requirement parameters. "Unclassified mapped" reads are those that mapped only to database sequences of species that are presumed to not be present in the sample by Emu's algorithm (likely due to low overall abundance).

```bash
#plot krona plots
sbatch ~/git_repos/Wrappers/unibz/run_krona.sh /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1/krona
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1/krona Control2 /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1/Control.txt
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1/krona Insecicide Insecicide.txt
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1/krona Microsap Microsap.txt

echo sample,group,color > metadata.csv
echo SHBB01881-1,Control,#1f77b4 >> metadata.csv
echo SHBB01881-2,Control,#1f77b4 >> metadata.csv
echo SHBB01881-3,Control,#1f77b4 >> metadata.csv
echo SHBB01884-3,Control,#1f77b4 >> metadata.csv
echo SHBB01887-2,Control,#1f77b4 >> metadata.csv
echo SHBB01890-1,Control,#1f77b4 >> metadata.csv
echo SHBB01895-1,Control,#1f77b4 >> metadata.csv
echo SHBB01900-1,Control,#1f77b4 >> metadata.csv
echo SHBB01888-1,Microsap,#2ca02c >> metadata.csv
echo SHBB01888-3,Microsap,#2ca02c >> metadata.csv
echo SHBB01888-4,Microsap,#2ca02c >> metadata.csv
echo SHBB01889-1,Microsap,#2ca02c >> metadata.csv
echo SHBB01891-1,Microsap,#2ca02c >> metadata.csv
echo SHBB01910-1,Microsap,#2ca02c >> metadata.csv
echo SHBB01914-1,Microsap,#2ca02c >> metadata.csv
echo SHBB01915-1,Microsap,#2ca02c >> metadata.csv
echo SHBB01882-1,Insecticide,#d62728 >> metadata.csv
echo SHBB01882-2,Insecticide,#d62728 >> metadata.csv
echo SHBB01882-4,Insecticide,#d62728 >> metadata.csv
echo SHBB01898-1,Insecticide,#d62728 >> metadata.csv
echo SHBB01899-1,Insecticide,#d62728 >> metadata.csv
echo SHBB01903-1,Insecticide,#d62728 >> metadata.csv
echo SHBB01909-1,Insecticide,#d62728 >> metadata.csv
echo SHBB01911-1,Insecticide,#d62728 >> metadata.csv

for taxa_level in species genus family order class phylum; do

#plot stacked barplots
module load apptainer/1.4.1-gcc-13.3.0-3coysxn
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_emu_abundance.py \
  -i /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 \
  -l "$taxa_level" \
  -n 9 \
  --output /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1/stacked-barplot_"$taxa_level".svg \
  --order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#plot heatmap plots
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_taxa_heatmap.py \
  -i /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 \
  -o /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 \
  -l "$taxa_level" \
  --log \
  --top_n 50 \
  --sample_order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#PCoA
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/emu_plot_pcoa.py \
  -i /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 \
  -o /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1 \
  -l "$taxa_level" \
  --metadata metadata.csv

done
```
***Repeat EMU with filtered reads***

```bash
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet-trim.csv
  Amplicon=16S
  Type=map-ont
#  Database=/data/users/theaven/db/emu/emu2026
#  Database=/data/users/theaven/db/emu/emu
#  Database=/data/users/theaven/db/emu/silva
  Database=/data/users/theaven/db/emu/silva-138.1
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/trim
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done

###################################################################################################

for Dir in $(ls -d /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/); do
  Task=emu
  InDir="$Dir"
  Sample_Sheet=/data/users/theaven/Ips_jam_project/20260323_sample_sheet-trim.csv
  Amplicon=ITS
  Type=map-ont
#  Database=/data/users/theaven/db/emu/silva
  Database=/data/users/theaven/db/emu/silva-138.1 #fails with no hits for SHBB01898-1
#  Database=/data/users/theaven/db/emu/unite-all
#  Database=/data/users/theaven/db/emu/unite-fungi
  Max_len=2000
  Min_len=0
  Abundance=0.0001
  OutDir=/data/users/theaven/Ips_jam_project/"$Task"/"$Amplicon"/"$(basename $Database)"/trim
  ExpectedOutput="$OutDir"/track.txt

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_EMU.sh "$InDir" "$Sample_Sheet" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" --abundance "$Abundance")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done
```
Plot for different runs:
```bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn
for InDir in $(ls -d /data/users/theaven/Ips_jam_project/emu/*/*/* | grep -v '/data/users/theaven/Ips_jam_project/emu/16S/emu2026/1' ); do
#plot krona plots
sbatch ~/git_repos/Wrappers/unibz/run_krona.sh "$InDir" "$InDir"/krona
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh "$InDir" "$InDir"/krona Control "$InDir"/Control.txt
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh "$InDir" "$InDir"/krona Insecicide "$InDir"/Insecicide.txt
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh "$InDir" "$InDir"/krona Microsap "$InDir"/Microsap.txt

for taxa_level in species genus family order class phylum; do

#plot stacked barplots
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_emu_abundance.py \
  -i "$InDir" \
  -l "$taxa_level" \
  -n 9 \
  --output "$InDir"/stacked-barplot_"$taxa_level".svg \
  --order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#plot heatmap plots
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_taxa_heatmap.py \
  -i "$InDir" \
  -o "$InDir" \
  -l "$taxa_level" \
  --log \
  --top_n 50 \
  --sample_order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#PCoA
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/emu_plot_pcoa.py \
  -i "$InDir" \
  -o "$InDir" \
  -l "$taxa_level" \
  --metadata /data/users/theaven/Ips_jam_project/metadata.csv

done
done
```

## NanoVI
NanoVI is a Nextflow DSL2 pipeline that performs taxonomic classification of full-length 16S ribosomal RNA gene sequences generated by Oxford Nanopore Technologies long-read sequencing. It uses variational inference to estimate species-level relative abundances from alignment likelihoods computed via CIGAR string analysis. 

From 10.64898/2026.03.07.710315 : NanoVI and Emu produced highly similar taxonomic profiles, successfully recovering all major taxa with appropriate relative abundances. In contrast, NanoCLUST and EPI2ME assigned a substantially larger fraction of reads to an "Other" category and failed to detect several key species.

NanoVI includes a helper script to build a reference database directly from the GTDB SSU FASTA file. The script assigns one taxid per species, deduplicates identical 16S sequences, and outputs the files required by the pipeline.

```bash
#NanoVI includes a helper script to build a reference database directly from the GTDB SSU FASTA file. The script assigns one taxid per species, deduplicates identical 16S sequences, and outputs the files required by the pipeline:
cd /home/clusterusers/theaven/db/GTDB
srun  -p bioagri  -c 1 --mem=4G --account=shame --pty bash 
module load apptainer/1.4.1-gcc-13.3.0-3  

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/build_gtdb_db.py \
    --ssu ./ssu_all_r226.fna.gz \
    --db-name db_gtdb_r226 \
    --output-dir . \
    --min-length 900

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/minimap2_2.30--h577a1d6_0 minimap2 -k 21 -d /home/clusterusers/theaven/db/GTDB/gtdb_index.mmi /home/clusterusers/theaven/db/GTDB/species_taxid.fasta

#Create the input samplesheet:              
echo sample,fastq > /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01881-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode01.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01881-2,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode02.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01881-3,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode03.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01882-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode04.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01882-2,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode05.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01882-4,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode06.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01884-3,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode07.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01887-2,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode08.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01888-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode09.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01888-3,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode10.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01888-4,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode11.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01889-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode12.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01890-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode13.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01891-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode14.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01895-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode15.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01898-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode16.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01899-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode17.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01900-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode18.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01903-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode19.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01909-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode20.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01910-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode21.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01911-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode22.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01914-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode23.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 
echo SHBB01915-1,/data/users/theaven/Ips_jam_project/raw_data/minion/16S/basecalls/barcode24.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv 

#Run nanovi:
for Input in $(ls -d /data/users/theaven/Ips_jam_project/nanovi/samplesheet.csv ); do
  Task=nanovi
  Amplicon=16S
  Type=map-ont
  Database=/home/clusterusers/theaven/db/GTDB
  Taxonomy=/home/clusterusers/theaven/db/GTDB/taxonomy.tsv
  Max_len=2000
  Min_len=500
  OutDir=/data/users/theaven/Ips_jam_project/nanovi/results/singularity
  ExpectedOutput="$OutDir"/consolidated_output/otu_table_abundance.tsv

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_nanovi.sh "$Input" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --taxonomy "$Taxonomy" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" )
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done
```
```bash
#plot krona plots
sbatch ~/git_repos/Wrappers/unibz/run_krona.sh /data/users/theaven/Ips_jam_project/nanovi/results/singularity/sample_outputs /data/users/theaven/Ips_jam_project/nanovi/results/singularity/krona
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh /data/users/theaven/Ips_jam_project/nanovi/results/singularity/sample_outputs /data/users/theaven/Ips_jam_project/nanovi/results/singularity/krona Control /data/users/theaven/Ips_jam_project/emu/16S/emu2026/1/Control.txt
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh /data/users/theaven/Ips_jam_project/nanovi/results/singularity/sample_outputs /data/users/theaven/Ips_jam_project/nanovi/results/singularity/krona Insecicide Insecicide.txt
sbatch ~/git_repos/Wrappers/unibz/run_krona_group.sh /data/users/theaven/Ips_jam_project/nanovi/results/singularity/sample_outputs /data/users/theaven/Ips_jam_project/nanovi/results/singularity/krona Microsap Microsap.txt

for taxa_level in species genus family order class phylum; do

#plot stacked barplots
module load apptainer/1.4.1-gcc-13.3.0-3coysxn
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_emu_abundance.py \
  -i /data/users/theaven/Ips_jam_project/nanovi/results/singularity/sample_outputs \
  -l "$taxa_level" \
  -n 9 \
  --output /data/users/theaven/Ips_jam_project/nanovi/results/singularity/stacked-barplot_"$taxa_level".svg \
  --order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#plot heatmap plots
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_taxa_heatmap.py \
  -i /data/users/theaven/Ips_jam_project/nanovi/results/singularity/sample_outputs \
  -o /data/users/theaven/Ips_jam_project/nanovi/results/singularity \
  -l "$taxa_level" \
  --log \
  --top_n 50 \
  --sample_order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#PCoA
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/emu_plot_pcoa.py \
  -i /data/users/theaven/Ips_jam_project/nanovi/results/singularity/sample_outputs \
  -o /data/users/theaven/Ips_jam_project/nanovi/results/singularity \
  -l "$taxa_level" \
  --metadata metadata.csv

done
```
***Repeat EMU with filtered reads***

```bash
#Create the input samplesheet:              
echo sample,fastq > /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01881-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode01.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01881-2,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode02.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01881-3,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode03.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01882-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode04.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01882-2,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode05.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01882-4,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode06.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01884-3,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode07.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01887-2,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode08.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01888-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode09.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01888-3,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode10.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01888-4,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode11.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01889-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode12.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01890-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode13.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01891-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode14.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01895-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode15.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01898-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode16.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01899-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode17.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01900-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode18.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01903-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode19.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01909-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode20.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01910-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode21.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01911-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode22.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01914-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode23.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 
echo SHBB01915-1,/data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/barcode24.trim.fastq >> /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv 

#Run nanovi:
for Input in $(ls -d /data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv ); do
  Task=nanovi
  Amplicon=16S
  Type=map-ont
  Database=/home/clusterusers/theaven/db/GTDB
  Taxonomy=/home/clusterusers/theaven/db/GTDB/taxonomy.tsv
  Max_len=2000
  Min_len=500
  OutDir=/data/users/theaven/Ips_jam_project/nanovi/results/singularity-trim
  ExpectedOutput="$OutDir"/consolidated_output/otu_table_abundance.tsv

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_nanovi.sh "$Input" "$OutDir" --Amplicon "$Amplicon" --database "$Database" --taxonomy "$Taxonomy" --type "$Type" --max_len "$Max_len" --min_len "$Min_len" )
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done
```

## UNIFRAC

UniFrac is a phylogeny-based distance metric, instead of just comparing abundance tables, it compares how evolutionarily related the organisms in each sample.

10.1111/1755-0998.13991

UniFrac depends on sequence-level features (ASVs/OTUs), not taxonomy
VSEARCH clustering at 97–99%
then build tree from representative OTUs

Nanopore reads
   ↓
ASV inference (DADA2 / VSEARCH / etc.)
   ↓
rep-seqs.qza
   ↓
MAFFT alignment
   ↓
IQ-TREE / FastTree
   ↓
UniFrac

Reads were dereplicated using usearch -fastx_uniques -minuniquesize 2 sizeout. ASVs were resolved using usearch -unoise3 -minsize 8.

Quality-filtered raw reads were de-replicated using USEARCH v11.0 (63) with the “-fastx_uniques” command and a minimum number of sequence observations of 2. ASVs were generated with the de-replicated sequences using the UNOISE3 algorithm (64) in USEARCH v11.0, with a minimum unique size of 10. OTUs clustered at 97% identity were generated from size-sorted sequences de-replicated with the “-cluster_otus” command in USEARCH v11.0.

The adaptor trimmed data was filtered using filtlong with the settings --min_length 3500  --min_mean_q 70 (v0.2.0 https://github.com/rrwick/Filtlong) and cutadapt3 (v2.7) with -m 3500  –M 6000.

For VSEARCH, sequences were processed using v.2.22.1  (Rognes et al., 2016). The CCS or UMI consensus sequences were  first dereplicated with the vsearch --derep_fulllength command  with the --sizeout and --relabel uniq options and then denoised  to resolve ASVs using first the vsearch -cluster_unoise command  with the --minsize option ranging from 2 to 8. Chimeras were finally detected and removed using the vsearch --uchime3_denovo  command.

Dorado (basecalling + adapter trimming)
→ Cutadapt (primer trimming)
→ Length filtering (e.g. 1200–1800 bp for 16S)
→ Optional: filtlong (moderate filtering)
→ EMU / NanoVI / wf-16s
#### Vsearch
10.1111/1755-0998.13991 use vsearch and usearch to genertate OTUs with unoise, others use usearch 10.1093/pnasnexus/pgae411 with unoise algorithms, they report good results with these denoising algorithms despite being designed for illumina data (DADA2 is bad) - these are with UMI nanopore data...
```bash
srun -p bioagri  -c 4 --mem 16G --pty bash
module load anaconda3
conda activate vsearch

#Dereplicate, denoise, de-chimera
#With minsize 2
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"

  vsearch -cluster_unoise "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --minsize 2 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"
done

for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"

  vsearch -cluster_unoise "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --minsize 2 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"
done

#With minsize 1
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"

  vsearch -cluster_unoise "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --minsize 1 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"
done

for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"

  vsearch -cluster_unoise "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --minsize 1 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"
done

for file in /data/users/theaven/Ips_jam_project/qc_data/minion/*/basecalls/CutAdapt/Filtlong/vsearch/*_1.log /data/users/theaven/Ips_jam_project/qc_data/minion/*/basecalls/CutAdapt/Filtlong/vsearch/*_2.log; do
OutFile=$(dirname $file)/minsize_$(basename $file | cut -d '_' -f2 | sed 's@.log@@g')_stats.tsv
if [[ ! -f "$OutFile" ]]; then
  echo -e "barcode\tsequences\tmin\tmax\tavg\tunique_sequences\tdiscarded_sequences\tpercentage_discarded\tclusters\tcluster_avg\tcluster_max\tsingletons\tpercent_seqs\tpercent_clusters\tmin\tmax\tavg\tchimeras\tnon_chimeras\ttotal_sequences" > "$OutFile"
fi
ID=$(basename "$file" | cut -d '_' -f1)
sequences=$(cat "$file" | grep -oP 'nt in \K[0-9]+' | head -1)
min=$(cat "$file" | grep -oP 'min \K[0-9]+' | head -1)
max=$(cat "$file" | grep -oP 'max \K[0-9]+' | head -1)
avg=$(cat "$file" | grep -oP 'avg \K[0-9]+' | head -1)
unique_sequences=$(cat "$file" | grep -oP '\d+(?= unique sequences)' | head -1) 
discarded_sequences=$(cat "$file" | grep -oP '\d+(?= sequences discarded)') 
discarded_sequences=${discarded_sequences:-0}
percentage_discarded=$(awk -v d="$discarded_sequences" -v u="$unique_sequences" '
BEGIN {
  if (u > 0) printf "%.2f", (d/u)*100;
  else print "NA"
}')
clusters=$(cat "$file" | grep -oP 'Clusters: \K[0-9]+')
cluster_avg=$(cat "$file" | grep -oP 'avg \K[0-9]+\.[0-9]+' | tail -1)
cluster_max=$(cat "$file" | grep -oP 'max \K[0-9]+' | tail -2 | head -1)
singletons=$(cat "$file" | grep -oP 'Singletons: \K[0-9]+')
percent_seqs=$(cat "$file" | grep -oP 'Singletons: [0-9]+, \K[0-9]+\.[0-9]+%')
percent_clusters=$(cat "$file" | grep -oP '\d+\.\d+%(?= of clusters)') 
min2=$(cat "$file" | grep -oP 'min \K[0-9]+' | tail -1)
max2=$(cat "$file" | grep -oP 'max \K[0-9]+' | tail -1)
avg2=$(cat "$file" | grep -oP 'avg \K[0-9]+' | tail -1)
chimeras=$(cat "$file" | grep -oP '\d+ \(\d+\.\d+%\)(?= chimeras)' | head -1) 
non_chimeras=$(cat "$file" | grep -oP 'chimeras, \K\d+ \(\d+\.\d+%\)' | head -1) 
total_sequences=$(cat "$file" | grep -oP 'in \K[0-9]+ total sequences' | grep -oP '[0-9]+')

echo -e "$ID\t$sequences\t$min\t$max\t$avg\t$unique_sequences\t$discarded_sequences\t$percentage_discarded\t$clusters\t$cluster_avg\t$cluster_max\t$singletons\t$percent_seqs\t$percent_clusters\t$min2\t$max2\t$avg2\t$chimeras\t$non_chimeras\t$total_sequences" >> "$OutFile"
done 
```
dataset had massive sequence-level noise
```bash
srun -p bioagri  -c 4 --mem 16G --pty bash
module load anaconda3
conda activate vsearch

#Dereplicate, denoise, de-chimera
#With minsize 2
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --sizeout --minuniquesize 2  \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.uc" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"
done

for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --sizeout --minuniquesize 2 \
    --relabel uniq  \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.uc" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"
done

#With minsize 1
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.uc" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"
done

for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.uc" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"
done

for file in /data/users/theaven/Ips_jam_project/qc_data/minion/*/basecalls/CutAdapt/Filtlong/vsearch/*cs.log; do
OutFile=$(dirname $file)/minsize_$(basename $file | cut -d '_' -f2 | sed 's@.log@@g')_stats.tsv
if [[ ! -f "$OutFile" ]]; then
  echo -e "barcode\tsequences\tmin\tmax\tavg\tunique_sequences\tdiscarded_sequences\tpercentage_discarded\tclusters\tcluster_avg\tcluster_max\tsingletons\tpercent_seqs\tpercent_clusters\tmin\tmax\tavg\tchimeras\tnon_chimeras\ttotal_sequences" > "$OutFile"
fi
ID=$(basename "$file" | cut -d '_' -f1)
sequences=$(cat "$file" | grep -oP 'nt in \K[0-9]+' | head -1)
min=$(cat "$file" | grep -oP 'min \K[0-9]+' | head -1)
max=$(cat "$file" | grep -oP 'max \K[0-9]+' | head -1)
avg=$(cat "$file" | grep -oP 'avg \K[0-9]+' | head -1)
unique_sequences=$(cat "$file" | grep -oP '\d+(?= unique sequences)' | head -1) 
discarded_sequences=$(cat "$file" | grep -oP '\d+(?= clusters discarded)') 
discarded_sequences=${discarded_sequences:-0}
percentage_discarded=$(awk -v d="$discarded_sequences" -v u="$unique_sequences" '
BEGIN {
  if (u > 0) printf "%.2f", (d/u)*100;
  else print "NA"
}')
clusters=$(cat "$file" | grep -oP 'Clusters: \K[0-9]+')
cluster_avg=$(cat "$file" | grep -oP 'avg \K[0-9]+\.[0-9]+' | tail -1)
cluster_max=$(cat "$file" | grep -oP 'max \K[0-9]+' | tail -2 | head -1)
singletons=$(cat "$file" | grep -oP 'Singletons: \K[0-9]+')
percent_seqs=$(cat "$file" | grep -oP 'Singletons: [0-9]+, \K[0-9]+\.[0-9]+%')
percent_clusters=$(cat "$file" | grep -oP '\d+\.\d+%(?= of clusters)') 
min2=$(cat "$file" | grep -oP 'min \K[0-9]+' | tail -1)
max2=$(cat "$file" | grep -oP 'max \K[0-9]+' | tail -1)
avg2=$(cat "$file" | grep -oP 'avg \K[0-9]+' | tail -1)
chimeras=$(cat "$file" | grep -oP '\d+ \(\d+\.\d+%\)(?= chimeras)' | head -1) 
non_chimeras=$(cat "$file" | grep -oP 'chimeras, \K\d+ \(\d+\.\d+%\)' | head -1) 
total_sequences=$(cat "$file" | grep -oP 'in \K[0-9]+ total sequences' | grep -oP '[0-9]+')

echo -e "$ID\t$sequences\t$min\t$max\t$avg\t$unique_sequences\t$discarded_sequences\t$percentage_discarded\t$clusters\t$cluster_avg\t$cluster_max\t$singletons\t$percent_seqs\t$percent_clusters\t$min2\t$max2\t$avg2\t$chimeras\t$non_chimeras\t$total_sequences" >> "$OutFile"
done 
```
Nanoclust
medaka/racon
```bash
#edit nextflow.config to contain conda.enabled = true

module load anaconda3
module load openjdk/17.0.11_9-none-none-2c62zhf
rm /home/clusterusers/theaven/tools/NanoCLUST/results/pipeline_info/execution_trace.txt
/home/clusterusers/theaven/tools/NanoCLUST/nextflow.1 run main.nf -profile test,conda
```
https://bugseq.com/free.

nanoASV
RAMBO
CONCOMPRA

- only takes one primer pair at once... this seems stupid design as the nanopore amplicon kit contains a primer mix and this will therefore be the situation for the vast majority of potential users.

primer-chop does not support degenerate bases
```bash

screen -S concompra
srun -p bioagri  -c 8 --mem 32G --pty bash
module load apptainer/1.4.1-gcc-13.3.0-3  
apptainer pull concompra.sif docker://willemstock/concompra:version0.0.2
mkdir /data/users/theaven/Ips_jam_project/concompra/ITS
cd /data/users/theaven/Ips_jam_project/concompra/ITS

#symlinked files not sufficient
cp /data/users/theaven/Ips_jam_project/raw_data/minion/ITS/basecalls/*.fastq /data/users/theaven/Ips_jam_project/concompra/ITS/.

echo TEMPLATE_DIR="/opt/CONCOMPRA/scripts" > /data/users/theaven/Ips_jam_project/concompra/ITS/directory_list.txt
echo PRIMER_SET="/data/users/theaven/Ips_jam_project/concompra/ITS/primer_set.fa" >> /data/users/theaven/Ips_jam_project/concompra/ITS/directory_list.txt
echo MIN=300 >> /data/users/theaven/Ips_jam_project/concompra/ITS/directory_list.txt
echo MAX=1500 >> /data/users/theaven/Ips_jam_project/concompra/ITS/directory_list.txt
echo MERGE_CONSENSUS=0.97 >> /data/users/theaven/Ips_jam_project/concompra/ITS/directory_list.txt
echo READS_CONSENSUS=120 >> /data/users/theaven/Ips_jam_project/concompra/ITS/directory_list.txt
echo THREADS=8 >> /data/users/theaven/Ips_jam_project/concompra/ITS/directory_list.txt

echo ">head" > /data/users/theaven/Ips_jam_project/concompra/ITS/primer_set.fa
echo TCCGTAGGTGAACCTGCGG >> /data/users/theaven/Ips_jam_project/concompra/ITS/primer_set.fa
echo ">tail" >> /data/users/theaven/Ips_jam_project/concompra/ITS/primer_set.fa
echo GCATATCAATAAGCGGAGGA >> /data/users/theaven/Ips_jam_project/concompra/ITS/primer_set.fa

#run the image, pointing to a local directory with has the (compressed) fastq files, the directory_list.txt (adjust the parameters but not the directories in this file) and the primer_set.fa (with the appropriate primer+anchor sequences) files
apptainer run --bind /data/users/theaven/Ips_jam_project/concompra/ITS:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/concompra.sif 
```
```bash
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime phylogeny

#demoising step:
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza

#phylogeny estimation
  qiime phylogeny align-to-tree-mafft-iqtree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-table table.qza \
  --i-phylogeny tree.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics

qiime diversity beta-phylogenetic \
  --i-table table.qza \
  --i-phylogeny tree.qza \
  --p-metric weighted_unifrac
```