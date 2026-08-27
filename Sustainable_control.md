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

#Fungal ITS sequences were downloaded from UNITE (DOI: https://doi.org/10.15156/BIO/3301229) [Accessed 1800202026; https://unite.ut.ee/repository.php], to this .fasta I will add the Ips sequences which have been predicted to be amplified by our ITS1/4 primmers
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

#62.4% of reads dropped, 8.18% because of no primer sequence in the reads, 0.29% too long, 51.73% tool short on average

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

#70.4% of reads dropped, 1.65% because of no primer sequence in the reads, 7.26% too long, 59.9% tool short on average

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
```bash
#plot stacked barplots
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

for file in $(ls /data/users/theaven/Ips_jam_project/wf-16s/out/*-trim/*/abundance_table_genus.tsv); do
  taxa_level=Family
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_epi2me_abundance2.py \
  -i "$file" \
  -n 9 \
  -o $(dirname "$file")/stacked-barplot_"$taxa_level".svg \
  --tax_level "$taxa_level" \
  --sample_order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01903-1 SHBB01909-1 SHBB01911-1
done

# -------------------------
# Plot
# -------------------------
fig, ax = plt.subplots(figsize=(16, 6))  # wider figure for more samples
pivot.plot(kind="bar", stacked=True, width=0.9, ax=ax)

plt.ylabel("Relative abundance")
plt.xlabel("Sample")
plt.xticks(rotation=45, ha="right")

# Legend top → bottom
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title=args.level,
          bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()

)

```

## EMU

Emu is a relative abundance estimator for 16S genomic sequences. The method is optimized for error-prone full-length reads, but can also be utilized for short-read data.

known to overinflate diversity, also reported to capture more taxa in mock dataset than alternatives.

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
  --order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#plot heatmap plots
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_taxa_heatmap.py \
  -i "$InDir" \
  -o "$InDir" \
  -l "$taxa_level" \
  --log \
  --top_n 50 \
  --sample_order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01903-1 SHBB01909-1 SHBB01911-1

#PCoA
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/emu_plot_pcoa.py \
  -i "$InDir" \
  -o "$InDir" \
  -l "$taxa_level" \
  --metadata /data/users/theaven/Ips_jam_project/metadata.csv

done
done

for InDir in $(ls -d /data/users/theaven/Ips_jam_project/emu/*/*/* | grep -v '/data/users/theaven/Ips_jam_project/emu/16S/emu2026/1' ); do
  for file in "$InDir"/*rel-abundance.tsv; do
    OutDir=$(dirname "$file")/4
mkdir "$OutDir"
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/propogate_taxa2.py -i "$file" -o "$OutDir"/$(basename "$file")
done
done

for InDir in $(ls -d /data/users/theaven/Ips_jam_project/emu/*/*/*/4 | grep -v '/data/users/theaven/Ips_jam_project/emu/16S/emu2026/1' ); do
  for taxa_level in species genus family order class phylum; do

#plot stacked barplots
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/plot_emu_abundance2.py \
  -i "$InDir" \
  -l "$taxa_level" \
  -n 13 \
  --output "$InDir"/stacked-barplot_"$taxa_level".svg \
  --order SHBB01881-1 SHBB01881-2 SHBB01881-3 SHBB01884-3 SHBB01887-2 SHBB01890-1 SHBB01895-1 SHBB01900-1 SHBB01882-1 SHBB01882-2 SHBB01882-4 SHBB01898-1 SHBB01899-1 SHBB01910-1 SHBB01914-1 SHBB01915-1 SHBB01888-1 SHBB01888-3 SHBB01888-4 SHBB01889-1 SHBB01891-1 SHBB01903-1 SHBB01909-1 SHBB01911-1
done
done
```
```bash
for InDir in /data/users/theaven/Ips_jam_project/emu/*/*/*/4; do
  cd "$InDir"
for f in *_rel-abundance.tsv; do
    sample=$(basename "$f" _rel-abundance.tsv)

    awk -v sample="$sample" '
    BEGIN { unmapped=0; total=0 }

    NR==1 { next }  # skip header

    {
        total += $10

        if ($1 == "unmapped") {
            unmapped += $10
        }
    }

    END {
        if (total > 0)
            printf "%s\t%.4f\t%.4f\t%.2f%%\n", sample, unmapped, total, (unmapped/total)*100
        else
            printf "%s\t0\t0\t0%%\n", sample
    }' "$f"

done > unmapped_percentage.tsv
done
```
***16S***

Stats and plotting in R - standardise with illumina workflow
```bash
cd /data/users/theaven/Ips_jam_project/emu/16S/silva/trim/4

for f in *_rel-abundance.tsv; do
  sample=${f/_rel-abundance.tsv/}

  awk -v s="$sample" '
  BEGIN {OFS="\t"}

  NR>1 {

    count = $10

    # keep only valid numeric counts
    if (count !~ /^[0-9.]+$/) next

    # build taxonomy string
    tax = $1 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9

    # filter junk taxa
    if (tax ~ /unidentified/ || tax ~ /;;;;/) next

    print tax, s, count
  }' "$f"

done \
| awk '
BEGIN {OFS="\t"}
{
  tax[$1]
  sam[$2]
  val[$1,$2]=$3
}

END {
  printf "taxon"
  for (s in sam) printf OFS s
  print ""

  for (t in tax) {
    printf t
    for (s in sam) {
      printf OFS (val[t,s]+0)
    }
    print ""
  }
}
' > EMU_otu_matrix_16s2.tsv
```
```R
long_otu <- read.table("EMU_otu_matrix_16s2.tsv", header=TRUE, sep="\t", row.names=1)
colnames(long_otu) <- gsub("\\.", "-", colnames(long_otu))
rownames(long_otu)[rownames(long_otu) == "mapped_unclassified;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root"] <- "mapped_unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified"
rownames(long_otu)[rownames(long_otu) == "mapped_filtered;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root"] <- "mapped_filtered;filtered;filtered;filtered;filtered;filtered;filtered;filtered"
rownames(long_otu)[rownames(long_otu) == "unmapped;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root"] <- "unmapped;unmapped;unmapped;unmapped;unmapped;unmapped;unmapped;unmapped"
new_names <- sub("^[^;]+;", "", rownames(long_otu))
long_otu$taxa <- new_names
long_otu_merged <- aggregate(
  . ~ taxa,
  data = long_otu,
  FUN = sum
)
rownames(long_otu_merged) <- long_otu_merged$taxa
long_otu_merged$taxa <- NULL
long_otu <- long_otu_merged

long_tax <- do.call(rbind, strsplit(rownames(long_otu), ";"))
colnames(long_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
long_tax <- long_tax[, 1:7]  # ensure correct columns only
tax_key <- apply(long_tax, 1, function(x) {
  paste(x, collapse = ";")
})
tax_key <- gsub("\\s+", "", tax_key)
rownames(long_tax) <- tax_key
```
Can be input to the illumina pipeline from here
```R
setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")
set.seed(1)

# Load metadata
meta <- read_tsv(
  "down_20260609/sample-metadata.tsv",
  comment = "",  
  show_col_types = FALSE
)
meta <- column_to_rownames(meta, var = "#SampleID")

# Create objects
TAX <- tax_table(as.matrix(long_tax))
OTU <- otu_table(as.matrix(long_otu), taxa_are_rows = TRUE)
SAM <- sample_data(meta)
ps <- phyloseq(OTU, SAM, TAX)

tax_table(ps) <- apply(tax_table(ps), 2, trimws)
tax_table(ps)[, "Genus"] <- gsub("^s__", "g__", tax_table(ps)[, "Genus"])

ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus <- subset_taxa(ps_genus, !Genus %in% c("unmapped", "filtered", "unclassified")) #remove unmapped

ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

df <- psmelt(ps_rel)

taxa_abund <- tapply(df$Abundance, df$Genus, sum)

top <- names(sort(taxa_abund, decreasing = TRUE))[1:13]

df$Genus <- as.character(df$Genus)
df$Genus[!df$Genus %in% top] <- "Other"

df$Genus <- factor(df$Genus, levels = c(top, "Other"))

auto_cols <- setNames(c(
  "#008000", "wheat3", "darkblue", "sienna3", "deeppink4", "#ff00ec","skyblue3", "orange" ,"#9467bf", "red3" ,"#71c837" ,"#000000", "#4dfad8"
), top[1:length(top)])

final_cols <- c(auto_cols, "Other" = "grey80")

df_sub <- subset(df, treatment != "blank")

df_sub$treatment <- as.character(df_sub$treatment)

df_sub$Group <- paste(df_sub$treatment, df_sub$Sample, sep = " ")
df_sub$Label <- paste(df_sub$treatment, df_sub$Sample)

df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

df_sub <- df_sub[order(df_sub$treatment, df_sub$Sample), ]

df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

ggplot(df_sub, aes(x = Label, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = final_cols) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
```
Alpha diversity:
```R
get_alpha <- function(ps_obj, meta_obj) {
  alpha <- estimate_richness(ps_obj,
                             measures = c("Shannon", "Simpson")) #chao is not appropriate for EMU proportional data and observed in sketchy
  rownames(alpha) <- gsub("\\.", "-", rownames(alpha))
  alpha$Sample <- rownames(alpha)
  meta_obj$Sample <- rownames(meta_obj)
  df <- merge(alpha, meta_obj, by = "Sample")
  return(df)
}

meta_f <- meta[meta$treatment != "blank", , drop = FALSE]  
ps_f <- prune_samples(rownames(meta_f), ps_genus)
meta_f <- meta_f[sample_names(ps_f), , drop = FALSE]
identical(sample_names(ps_f), rownames(meta_f))
alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = treatment, y = Shannon, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = treatment),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  scale_fill_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")

#######

by(alpha_df$Shannon, alpha_df$treatment, shapiro.test)
#Control - p-value = 0.6153 - normal
#Insecticide - p-value = 0.6341 - normal
#Microsap - p-value = 0.1651 - normal

pairwise.t.test(alpha_df$Shannon,
                alpha_df$treatment,
                p.adjust.method = "BH",
                pool.sd = FALSE)
#            Control Insecticide
#Insecticide 0.38    -          
#Microsap    0.51    0.31       

pairwise.wilcox.test(alpha_df$Shannon,
                     alpha_df$treatment,
                     p.adjust.method = "BH")
#No significant differences in Shannon diversity between any pair of treatments
#            Control Insecticide
#Insecticide 0.57    -          
#Microsap    0.65    0.48 
```
beta diversity:

```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    max.overlaps = Inf,
    force = 2
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

#######

betadisper_res <- betadisper(bray, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.00869 0.0043441 0.2625    999  0.773
#Residuals 21 0.34756 0.0165503 
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(bray ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.5000 0.05983 0.6682  0.888
#Residual 21   7.8558 0.94017              
#Total    23   8.3557 1.00000  
#No significant difference in community composition between treatments, treatment explains only ~5.9% of variation in community composition     

#######

ps_simple <- phyloseq::phyloseq(
  phyloseq::otu_table(ps_f),
  phyloseq::sample_data(ps_f)
)
ps_pa <- transform_sample_counts(ps_simple, function(x) as.numeric(x > 0))
jaccard <- phyloseq::distance(ps_pa, method = "jaccard")
ord <- ordinate(ps_f, method = "PCoA", distance = jaccard)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")


  ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

#######


betadisper_res <- betadisper(jaccard, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.004472 0.0022358 0.2319    999  0.775
#Residuals 21 0.202460 0.0096410
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(jaccard ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.7645 0.10106 1.1805    0.2
#Residual 21   6.7997 0.89894              
#Total    23   7.5641 1.00000 
#No significant difference in community composition between treatments, treatment explains ~10.1% of variation     
```
No differences

Presence/abscence - bubble plot:
```R
#build matrix for Genus level
ps_f_genus <- tax_glom(ps_f, taxrank = "Genus")
ps_f_genus <- filter_taxa(ps_f_genus, function(x) sum(x) > 0, TRUE)

mat <- as(otu_table(ps_f_genus), "matrix")
if (taxa_are_rows(ps_f_genus)) {
  mat <- t(mat)
}

pa_mat <- (mat > 0) * 1

meta_mat <- meta_f[match(rownames(pa_mat), rownames(meta_f)), , drop = FALSE]

group <- meta_mat$treatment
groups <- unique(group)
stopifnot(all(rownames(pa_mat) == rownames(meta_mat)))

#Pairwise Fisher tests (ALL pairs)
pair_list <- combn(groups, 2, simplify = FALSE)

pairwise_p <- lapply(pair_list, function(grp) {
  g1 <- grp[1]
  g2 <- grp[2]
  idx <- group %in% c(g1, g2)
  apply(pa_mat[idx, , drop = FALSE], 2, function(x) {
    tab <- table(x, group[idx])
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA)
    if (all(tab == 0)) return(NA)
    fisher.test(tab)$p.value
  })
})

names(pairwise_p) <- sapply(pair_list, paste, collapse = "_vs_")

pairwise_p_adj <- lapply(pairwise_p, function(p) {
  p.adjust(p, method = "BH")
})

#Extract significant taxa (ANY comparison)
sig_taxa <- unique(unlist(lapply(pairwise_p_adj, function(p) {
  names(p)[which(p < 0.05 & !is.na(p))]
})))

prev_list <- lapply(groups, function(g) {
  colMeans(pa_mat[group == g, sig_taxa, drop = FALSE])
})
names(prev_list) <- groups

tax <- as.data.frame(tax_table(ps_f_genus))
tax_labels <- tax$Genus
names(tax_labels) <- rownames(tax)

df <- data.frame(
  Taxon = sig_taxa,
  Label = tax_labels[sig_taxa]
)

for (g in groups) {
  df[[paste0(g, "_prev")]] <- prev_list[[g]]
}

#Long format for plotting
df_long <- pivot_longer(
  df,
  cols = ends_with("_prev"),
  names_to = "Group",
  values_to = "Prevalence"
)

df_long$Group <- gsub("_prev", "", df_long$Group)

df_long$Label[is.na(df_long$Label)] <- df_long$Taxon

ggplot(df_long, aes(
  x = Group,
  y = Label,
  size = Prevalence,
  color = Group
)) +
  geom_point(alpha = 0.85) +
  scale_size(range = c(2, 10)) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Genus",
    size = "Prevalence",
    color = "Group",
    title = "Significant taxa (pairwise Fisher tests, FDR < 0.05)"
  )
```
No significant taxa - hard with only 8 reps

Presence/abscence + adundance - heatmap:
```R
tax_tab <- as.data.frame(tax_table(ps_f_genus))

tax_names <- tax_tab$Genus
tax_names[is.na(tax_names) | tax_names == ""] <- "Unknown"
tax_names <- make.unique(tax_names)

mat_asv <- mat
colnames(mat_asv) <- tax_names

mat_genus_rel <- sweep(mat_asv, 1, rowSums(mat_asv), "/")
mat_genus_rel_log <- log10(mat_genus_rel + 1e-6)

#Order samples
meta_mat <- meta[rownames(mat_genus_rel_log), , drop = FALSE ]

ord <- order(meta_mat$treatment)

mat_ordered <- mat_genus_rel_log[ord, ]
meta_ordered <- meta_mat[ord, , drop = FALSE  ]

# --- Define colors: 0 = white, then blue → red ---
# Avoid including 0 in gradient
nonzero_vals <- mat_ordered[mat_ordered > 0]

# --- Row annotations ---
annotation_row <- data.frame(
  Treatment = meta_ordered$treatment
)

rownames(annotation_row) <- rownames(mat_ordered)

#Optional: gaps
gaps <- cumsum(table(meta_ordered$treatment))

pheatmap(
  mat_ordered,
  color = colorRampPalette(c("white", "blue", "red"))(100),
  breaks = seq(
    min(mat_ordered, na.rm = TRUE),
    max(mat_ordered, na.rm = TRUE),
    length.out = 101
  ),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  gaps_row = gaps,
  border_color = "grey90",
  fontsize_col = 5
)
```

***ITS***

Stats and plotting in R - standardise with illumina workflow
```bash
cd /data/users/theaven/Ips_jam_project/emu/ITS/unite-all/trim/4

for f in *_rel-abundance.tsv; do
  sample=${f/_rel-abundance.tsv/}

  awk -v s="$sample" '
  BEGIN {OFS="\t"}

  NR>1 {

    count = $10

    # keep only valid numeric counts
    if (count !~ /^[0-9.]+$/) next

    # build taxonomy string
    tax = $1 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9

    # filter junk taxa
    if (tax ~ /unidentified/ || tax ~ /;;;;/) next

    print tax, s, count
  }' "$f"

done \
| awk '
BEGIN {OFS="\t"}
{
  tax[$1]
  sam[$2]
  val[$1,$2]=$3
}

END {
  printf "taxon"
  for (s in sam) printf OFS s
  print ""

  for (t in tax) {
    printf t
    for (s in sam) {
      printf OFS (val[t,s]+0)
    }
    print ""
  }
}
' > EMU_otu_matrix_its2.tsv
```
```R
long_otu <- read.table("EMU_otu_matrix_its2.tsv", header=TRUE, sep="\t", row.names=1)
colnames(long_otu) <- gsub("\\.", "-", colnames(long_otu))
rownames(long_otu)[rownames(long_otu) == "mapped_unclassified;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root"] <- "mapped_unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified"
rownames(long_otu)[rownames(long_otu) == "mapped_filtered;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root"] <- "mapped_filtered;filtered;filtered;filtered;filtered;filtered;filtered;filtered"
rownames(long_otu)[rownames(long_otu) == "unmapped;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root"] <- "unmapped;unmapped;unmapped;unmapped;unmapped;unmapped;unmapped;unmapped"
rownames(long_otu)[rownames(long_otu) == "122;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root"] <- "unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root;unclassified_root2"
new_names <- sub("^[^;]+;", "", rownames(long_otu))
long_otu$taxa <- new_names
long_otu_merged <- aggregate(
  . ~ taxa,
  data = long_otu,
  FUN = sum
)
rownames(long_otu_merged) <- long_otu_merged$taxa
long_otu_merged$taxa <- NULL
long_otu <- long_otu_merged

long_tax <- do.call(rbind, strsplit(rownames(long_otu), ";"))
colnames(long_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
long_tax <- long_tax[, 1:7]  # ensure correct columns only
tax_key <- apply(long_tax, 1, function(x) {
  paste(x, collapse = ";")
})
tax_key <- gsub("\\s+", "", tax_key)
rownames(long_tax) <- tax_key
```
Can be input to the illumina pipeline from here
```R
setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")
set.seed(1)

# Load metadata
meta <- read_tsv(
  "down_20260609/sample-metadata.tsv",
  comment = "",  
  show_col_types = FALSE
)
meta <- column_to_rownames(meta, var = "#SampleID")

# Create objects
TAX <- tax_table(as.matrix(long_tax))
OTU <- otu_table(as.matrix(long_otu), taxa_are_rows = TRUE)
SAM <- sample_data(meta)
ps <- phyloseq(OTU, SAM, TAX)

tax_table(ps) <- apply(tax_table(ps), 2, trimws)
tax_table(ps)[, "Genus"] <- gsub("^s__", "g__", tax_table(ps)[, "Genus"])

ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus <- subset_taxa(ps_genus, !Genus %in% c("unmapped", "filtered", "unclassified")) #remove unmapped

ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

df <- psmelt(ps_rel)

taxa_abund <- tapply(df$Abundance, df$Genus, sum)

top <- names(sort(taxa_abund, decreasing = TRUE))[1:13]

df$Genus <- as.character(df$Genus)
df$Genus[!df$Genus %in% top] <- "Other"

df$Genus <- factor(df$Genus, levels = c(top, "Other"))

auto_cols <- setNames(c(
  "darkblue", "deeppink4", "orange","#00ffff", "sienna3", "wheat3","#ccffaa" ,"#00ffcc" ,"yellow", "#6ea02c", "#008000", "red3", "#ff00ec"
), top[1:length(top)])

final_cols <- c(auto_cols, "Other" = "grey80")

df_sub <- subset(df, treatment != "blank")

df_sub$treatment <- as.character(df_sub$treatment)

df_sub$Group <- paste(df_sub$treatment, df_sub$Sample, sep = " ")
df_sub$Label <- paste(df_sub$treatment, df_sub$Sample)

df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

df_sub <- df_sub[order(df_sub$treatment, df_sub$Sample), ]

df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

ggplot(df_sub, aes(x = Label, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = final_cols) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
```
Alpha diversity:
```R
get_alpha <- function(ps_obj, meta_obj) {
  alpha <- estimate_richness(ps_obj,
                             measures = c("Shannon", "Simpson")) #chao is not appropriate for EMU proportional data and observed in sketchy
  rownames(alpha) <- gsub("\\.", "-", rownames(alpha))
  alpha$Sample <- rownames(alpha)
  meta_obj$Sample <- rownames(meta_obj)
  df <- merge(alpha, meta_obj, by = "Sample")
  return(df)
}

meta_f <- meta[meta$treatment != "blank", , drop = FALSE]  
ps_f <- prune_samples(rownames(meta_f), ps_genus)
meta_f <- meta_f[sample_names(ps_f), , drop = FALSE]
identical(sample_names(ps_f), rownames(meta_f))
alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = treatment, y = Shannon, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = treatment),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  scale_fill_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")

#######

by(alpha_df$Shannon, alpha_df$treatment, shapiro.test)
#Control - p-value = 0.01199 - not normal
#Insecticide - p-value = 0.5794 - normal
#Microsap - p-value = 0.9767 - normal

pairwise.wilcox.test(alpha_df$Shannon,
                     alpha_df$treatment,
                     p.adjust.method = "BH")
#No significant differences in Shannon diversity between any pair of treatments
#            Control Insecticide
#Insecticide 0.35    -          
#Microsap    0.25    0.96 
```
beta diversity:

```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    max.overlaps = Inf,
    force = 2
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

#######

betadisper_res <- betadisper(bray, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.09413 0.047067 1.5849    999  0.228
#Residuals 21 0.62365 0.029698 
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(bray ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.6562 0.09005 1.0391  0.384
#Residual 21   6.6307 0.90995              
#Total    23   7.2869 1.00000   
#No significant difference in community composition between treatments, treatment explains only ~9.0% of variation in community composition     

#######

ps_simple <- phyloseq::phyloseq(
  phyloseq::otu_table(ps_f),
  phyloseq::sample_data(ps_f)
)
ps_pa <- transform_sample_counts(ps_simple, function(x) as.numeric(x > 0))
jaccard <- phyloseq::distance(ps_pa, method = "jaccard")
ord <- ordinate(ps_f, method = "PCoA", distance = jaccard)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")


  ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

#######


betadisper_res <- betadisper(jaccard, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.003626 0.0018129 0.2333    999  0.794
#Residuals 21 0.163153 0.0077692 
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(jaccard ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.5180 0.08328 0.9539  0.558
#Residual 21   5.7016 0.91672              
#Total    23   6.2195 1.00000 
#No significant difference in community composition between treatments, treatment explains ~8.3% of variation     
```
No differences

Presence/abscence - bubble plot:
```R
#build matrix for Genus level
ps_f_genus <- tax_glom(ps_f, taxrank = "Genus")
ps_f_genus <- filter_taxa(ps_f_genus, function(x) sum(x) > 0, TRUE)

mat <- as(otu_table(ps_f_genus), "matrix")
if (taxa_are_rows(ps_f_genus)) {
  mat <- t(mat)
}

pa_mat <- (mat > 0) * 1

meta_mat <- meta_f[match(rownames(pa_mat), rownames(meta_f)), , drop = FALSE]

group <- meta_mat$treatment
groups <- unique(group)
stopifnot(all(rownames(pa_mat) == rownames(meta_mat)))

#Pairwise Fisher tests (ALL pairs)
pair_list <- combn(groups, 2, simplify = FALSE)

pairwise_p <- lapply(pair_list, function(grp) {
  g1 <- grp[1]
  g2 <- grp[2]
  idx <- group %in% c(g1, g2)
  apply(pa_mat[idx, , drop = FALSE], 2, function(x) {
    tab <- table(x, group[idx])
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA)
    if (all(tab == 0)) return(NA)
    fisher.test(tab)$p.value
  })
})

names(pairwise_p) <- sapply(pair_list, paste, collapse = "_vs_")

pairwise_p_adj <- lapply(pairwise_p, function(p) {
  p.adjust(p, method = "BH")
})

#Extract significant taxa (ANY comparison)
sig_taxa <- unique(unlist(lapply(pairwise_p_adj, function(p) {
  names(p)[which(p < 0.05 & !is.na(p))]
})))

prev_list <- lapply(groups, function(g) {
  colMeans(pa_mat[group == g, sig_taxa, drop = FALSE])
})
names(prev_list) <- groups

tax <- as.data.frame(tax_table(ps_f_genus))
tax_labels <- tax$Genus
names(tax_labels) <- rownames(tax)

df <- data.frame(
  Taxon = sig_taxa,
  Label = tax_labels[sig_taxa]
)

for (g in groups) {
  df[[paste0(g, "_prev")]] <- prev_list[[g]]
}

#Long format for plotting
df_long <- pivot_longer(
  df,
  cols = ends_with("_prev"),
  names_to = "Group",
  values_to = "Prevalence"
)

df_long$Group <- gsub("_prev", "", df_long$Group)

df_long$Label[is.na(df_long$Label)] <- df_long$Taxon

ggplot(df_long, aes(
  x = Group,
  y = Label,
  size = Prevalence,
  color = Group
)) +
  geom_point(alpha = 0.85) +
  scale_size(range = c(2, 10)) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Genus",
    size = "Prevalence",
    color = "Group",
    title = "Significant taxa (pairwise Fisher tests, FDR < 0.05)"
  )
```
No significant taxa - hard with only 8 reps

Presence/abscence + adundance - heatmap:
```R
tax_tab <- as.data.frame(tax_table(ps_f_genus))

tax_names <- tax_tab$Genus
tax_names[is.na(tax_names) | tax_names == ""] <- "Unknown"
tax_names <- make.unique(tax_names)

mat_asv <- mat
colnames(mat_asv) <- tax_names

mat_genus_rel <- sweep(mat_asv, 1, rowSums(mat_asv), "/")
mat_genus_rel_log <- log10(mat_genus_rel + 1e-6)

#Order samples
meta_mat <- meta[rownames(mat_genus_rel_log), , drop = FALSE ]

ord <- order(meta_mat$treatment)

mat_ordered <- mat_genus_rel_log[ord, ]
meta_ordered <- meta_mat[ord, , drop = FALSE  ]

# --- Define colors: 0 = white, then blue → red ---
# Avoid including 0 in gradient
nonzero_vals <- mat_ordered[mat_ordered > 0]

# --- Row annotations ---
annotation_row <- data.frame(
  Treatment = meta_ordered$treatment
)

rownames(annotation_row) <- rownames(mat_ordered)

#Optional: gaps
gaps <- cumsum(table(meta_ordered$treatment))

pheatmap(
  mat_ordered,
  color = colorRampPalette(c("white", "blue", "red"))(100),
  breaks = seq(
    min(mat_ordered, na.rm = TRUE),
    max(mat_ordered, na.rm = TRUE),
    length.out = 101
  ),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  gaps_row = gaps,
  border_color = "grey90",
  fontsize_col = 7
)
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
***Repeat nanovi with filtered reads***

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
#With minsize 1
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
    --minsize 1 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.nonchimera.fa" \
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
    --minsize 1 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1.log"
done

#With minsize 2
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
    --minsize 2 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.nonchimera.fa" \
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
    --minsize 2 --id 0.97 \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2.nonchimera.fa" \
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

Try clustering with course grouping by similarity alone - no denoising.
```bash
srun -p bioagri  -c 4 --mem 16G --pty bash
module load anaconda3
conda activate vsearch

#Dereplicate, de-chimera
#With minsize 1
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.uc" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"
done

for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --sizeout \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.uc" \
    --threads 4 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-1cs.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_1cs.log"
done

#With minsize 2
for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --sizeout --minuniquesize 2  \
    --relabel uniq \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.uc" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.nonchimera.fa" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"
done

for Reads in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.fastq.gz; do
  OutDir="$(dirname "$Reads")/vsearch"
  mkdir -p "$OutDir"
  seqtk seq -A "$Reads" > tmp.fasta
  vsearch --derep_fulllength tmp.fasta \
    --output "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --sizeout --minuniquesize 2 \
    --relabel uniq  \
    --threads 1 2>&1 | tee "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --cluster_size "$OutDir/$(basename "$Reads" .fastq.gz).unique-2cs.fa" \
    --id 0.97 --strand both --sizein --sizeout \
    --centroids "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --uc "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.uc" \
    --threads 1 2>&1 | tee -a "$OutDir/$(basename "$Reads" .trim.filtlong.fastq.gz)_2cs.log"

  vsearch --uchime3_denovo "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.fa" \
    --nonchimeras "$OutDir/$(basename "$Reads" .fastq.gz).centroids-2cs.nonchimera.fa" \
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
Still high levels of noise - many OTUs with minsize 1, with minsize 2 >95% of reads are dropped due to high error rate of nanopore reads.

#### Kraken

Kraken has been run within the EPI2ME pipeline - this missassigned arthopoda reads based on the results of the other classigfication tools. Will try running seperetedly with different databases.

```bash
for Reads in $(ls /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.trim.filtlong.fastq.gz | grep -v 'unclassified\|code24'); do
Task=kraken
Database=/data/kraken2/pluspf
OutDir="$(dirname "$Reads")/kraken2.17.1/$(basename  "$Database")"
Samplesheet=/data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv
OutPrefix=$(awk -F',' -v f="$(basename "$Reads" .trim.filtlong.fastq.gz)" '$2 ~ "/ITS/" && $2 ~ f {print $1}' "$Samplesheet")
ID="$OutPrefix"_"$(basename  "$Database")"
ExpectedOutput="$OutDir"/"$OutPrefix"_report.txt

  Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  while [ "$Jobs" -gt 1 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  done

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable /data/users/theaven/run_kraken2.sh "$Reads" "$Database" "$OutDir" "$OutPrefix" )
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi

done

for Reads in $(ls /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/Filtlong/*.trim.filtlong.fastq.gz | grep -v 'unclassified\|code24'); do
Task=kraken
Database=/data/kraken2/NCBI/k2_NCBI_reference_20251007
OutDir="$(dirname "$Reads")/kraken2.17.1/$(basename  "$Database")"
Samplesheet=/data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv
OutPrefix=$(awk -F',' -v f="$(basename "$Reads" .trim.filtlong.fastq.gz)" '$2 ~ "/ITS/" && $2 ~ f {print $1}' "$Samplesheet")
ID="$OutPrefix"_"$(basename  "$Database")"
ExpectedOutput="$OutDir"/"$OutPrefix"_report.txt

  Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  while [ "$Jobs" -gt 1 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  done

  if [ ! -s "$ExpectedOutput" ]; then
    mkdir -p $OutDir
    ID=$(echo "$OutDir" | rev | cut -d '/' -f1 | rev)

    jobid=$(sbatch --job-name="$Task" --parsable /data/users/theaven/run_kraken2.sh "$Reads" "$Database" "$OutDir" "$OutPrefix" )
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi

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

### BLASTN <a name="12"></a>

BLAST is alignment based, limited by 'best hit' interpretation, and prone to misidentifying short or conserved seqeunces. For eDNA/metabarcoding, BLAST is often too literal — it finds the closest sequence, even if it’s wrong. However, the NCBI nt database is far larger than any dedicated database, including SILVA.

BLASTN vs NCBI nt database:
```bash
Task=blast
  Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  while [ "$Jobs" -gt 3 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  done

  #Run BLAST with many hits
for ASV in $(find /data/users/theaven/Ips_jam_project/qc_data/minion/16S/basecalls/CutAdapt/ -name '*.trim.fastq' -type f | grep -v 'unclassified'); do
  Task=blast
  Database=/data/blobtoolkit/nt/nt
  Max_target=10
  Samplesheet=/data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv
  OutPrefix=$(awk -F',' -v f="$ASV" '$2==f {print $1}' "$Samplesheet")
  seqkit fq2fa "$ASV" > "$(dirname $ASV)"/"$OutPrefix".fasta
  OutDir="$(dirname $ASV)"/"$Task"
  echo "$OutPrefix"
  mkdir -p $OutDir
  ExpectedOutput="$OutDir"/${OutPrefix}.vs."$(basename $Database)".mts"$Max_target".hsp1.1e25.megablast.out

  Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  while [ "$Jobs" -gt 4 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  done

  if [ ! -s "$ExpectedOutput" ]; then
    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_blastn2.sh "$(dirname $ASV)"/"$OutPrefix".fasta "$Database" "$OutDir" "$OutPrefix" "$Max_target")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done
```
Get LCA of BLAST hits for each read and collate into abundance table:
```bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

for file in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/*.vs.nt.mts10.hsp1.1e25.megablast.out; do
  mkdir $(dirname "$file")/LCA
awk 'BEGIN{OFS="\t"} {
print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$3,$2,$14
}' "$file" | sed '1d' > $(dirname "$file")/LCA/$(basename "$file")

cuttoff_rank=family 
apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/ncbi_blast_2_lca.py \
  --input $(dirname "$file")/LCA/$(basename "$file") \
  --nodes /data/users/theaven/db/blastn/taxonomy_db/nodes.dmp \
  --names /data/users/theaven/db/blastn/taxonomy_db/names.dmp \
  --outdir $(dirname "$file")/LCA \
  --min-rank "$cuttoff_rank"

sed -i 's/ /_/g' $(dirname "$file")/LCA/$(basename "$file" | sed 's@.out@@g').min-"$cuttoff_rank".lca.tsv

awk 'NR>1 {sum += $8} END {print sum}' $(dirname "$file")/LCA/$(basename "$file" | sed 's@.out@@g').min-"$cuttoff_rank".lca.tsv
done
```
```bash
module load anaconda3
conda activate taxonkit
for file in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/SHBB01898-1.vs.nt.mts10.hsp1.1e25.megablast.out; do

  outdir=$(dirname "$file")/lca
  mkdir -p "$outdir"

  base=$(basename "$file" .out)
  tmp="$outdir/${base}.tsv"

  # 1. top 10 hits per read (by bitscore)
  sort -k1,1 -k3,3nr "$file" | \
  awk '{
    if(count[$1] < 10) {
      print;
      count[$1]++;
    }
  }' | \
  cut -f1,2 | \
  sed '/^\s*$/d' > "$tmp"

  # 2. LCA
  taxonkit lca "$tmp" > "$outdir/${base}_lca.tsv"

done

cut -f1,2 /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/SHBB01898-1.vs.nt.mts10.hsp1.1e25.megablast.out > /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/lca/SHBB01898-1.vs.nt.mts10.hsp1.1e25.megablast.tsv
taxonkit lca /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/lca/SHBB01898-1.vs.nt.mts10.hsp1.1e25.megablast.tsv
```
```bash
module load anaconda3
conda activate megan

export _JAVA_OPTIONS="-Djava.io.tmpdir=$HOME/tmp"
mkdir -p $HOME/tmp

for file in /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/*.vs.nt.mts10.hsp1.1e25.megablast.out; do
  mkdir $(dirname "$file")/megan
awk 'BEGIN{OFS="\t"} {
print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$3,$2,$14
}' "$file" | sed '1d' > $(dirname "$file")/megan/$(basename "$file")

blast2lca \
  -i $(dirname "$file")/megan/$(basename "$file") \
  -f BlastTab -m BlastN \
  -a2t /data/users/theaven/db/ncbi/nucl_gb.accession2taxid \
  -top 10 \
  -ms 50 \
  -me 0.01 \
  -o $(dirname "$file")/megan/$(basename "$file")_lca.txt
done

#  -mdb /data/users/theaven/db/megan/megan-nr-r2.mdb \






awk 'NR>1 {sum += $8} END {print sum}' /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/LCA/id_counts.tsv

awk 'NR>1 {sum += $8} END {print sum}' temp.tsv

awk -F'\t' 'NR > 1 { 
    key = $1; 
    for (i=2; i<=7; i++) key = key "\t" $i; 
    counts[key]++; 
    total++ 
} 
END { 
    for (k in counts) { 
        printf "%s\t%d\t%.8f\n", k, counts[k], counts[k]/total 
    } 
}' /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/LCA/SHBB01903-1.vs.nt.mts10.hsp1.1e25.megablast.min-family.LCA_per_read.debug.tsv | sort -t$'\t' -k8,8nr
```
```python
from collections import defaultdict

infile = "/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/LCA/SHBB01903-1.vs.nt.mts10.hsp1.1e25.megablast.min-family.LCA_per_read.debug.tsv"

counts = defaultdict(int)
total = 0

with open(infile, "r") as f:
    header = next(f)  # skip header
    for line in f:
        if not line.strip():
            continue
        total += 1
        cols = line.rstrip("\n").split("\t")
        key = tuple(cols[:7])
        counts[key] += 1

with open("id_counts.tsv", "w") as out:
    for key, c in counts.items():
        proportion = c / total if total else 0
        out.write("\t".join(map(str, key)) + f"\t{c}\t{proportion}\n")
````
```python
import pandas as pd
from collections import defaultdict

# -----------------------------
# 1. Load BLAST file
# -----------------------------
cols = [
    "read_id", "subject", "pid", "aln_len",
    "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send",
    "evalue", "bit_score",
    "taxid", "taxon"
]

df = pd.read_csv(
    "/data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/blast/megan/SHBB01898-1.vs.nt.mts10.hsp1.1e25.megablast.out",
    sep="\t",
    names=cols
)

df["taxid"] = df["taxid"].astype(str).str.split(";").str[0].astype(int)

# -----------------------------
# 2. Load NCBI taxonomy (nodes.dmp)
# -----------------------------
# Format: taxid -> parent_taxid

def load_nodes(nodes_file):
    parent = {}
    with open(nodes_file, "r") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            parent_taxid = int(parts[1])
            parent[taxid] = parent_taxid
    return parent

parent_map = load_nodes("/data/users/theaven/db/blastn/taxonomy_db/nodes.dmp")

# -----------------------------
# 3. Build ancestry path
# -----------------------------
def get_ancestors(taxid, parent_map):
    path = set()
    while taxid in parent_map and taxid != parent_map[taxid]:
        path.add(taxid)
        taxid = parent_map[taxid]
    return path

# -----------------------------
# 4. MEGAN-style hit filtering
#    (keep hits within 90% of best bit score)
# -----------------------------
def filter_hits(group):
    best = group["bit_score"].max()
    threshold = best * 0.90
    return group[group["bit_score"] >= threshold]

filtered = df.groupby("read_id", group_keys=False).apply(filter_hits)

# -----------------------------
# 5. Lowest Common Ancestor (LCA)
# -----------------------------
def lca_taxids(taxids, parent_map):
    paths = []
    for t in taxids:
        paths.append(get_ancestors(t, parent_map))
    if not paths:
        return None
    common = set.intersection(*paths)
    return max(common) if common else 1  # fallback: root

def assign_lca(group):
    taxids = group["taxid"].unique()
    return lca_taxids(taxids, parent_map)

lca_series = filtered.groupby("read_id").apply(assign_lca)

# -----------------------------
# 6. Load names.dmp (taxid -> name)
# -----------------------------
def load_names(names_file):
    names = {}
    with open(names_file, "r") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            name = parts[1]
            unique_name = parts[2]
            name_class = parts[3]
            if "scientific name" in name_class:
                names[taxid] = name
    return names

names_map = load_names("/data/users/theaven/db/blastn/taxonomy_db/names.dmp")

# -----------------------------
# 7. Convert LCA taxid → name
# -----------------------------
lca_names = lca_series.map(lambda x: names_map.get(x, f"taxid_{x}"))

# -----------------------------
# 8. Abundance table
# -----------------------------
abundance = lca_names.value_counts().reset_index()
abundance.columns = ["taxon", "read_count"]

abundance["relative_abundance"] = (
    abundance["read_count"] / abundance["read_count"].sum()
)

# -----------------------------
# 9. Save results
# -----------------------------
abundance.to_csv("lca_abundance.tsv", sep="\t", index=False)

print(abundance.head(20))
```

### ITSx

Luciano suggests splitting the long ITS reads into constituent ITS1 and ITS2 regions and seeing if these can be classified better.

```bash
for file in $(find /data/users/theaven/Ips_jam_project/qc_data/minion/ITS/basecalls/CutAdapt/ -name '*.trim.fastq' -type f | grep -v 'unclassified'); do
  Task=ITSx
  Samplesheet=/data/users/theaven/Ips_jam_project/nanovi/samplesheet-trim.csv
  OutPrefix=$(awk -F',' -v f="$file" '$2==f {print $1}' "$Samplesheet")
  OutDir="$(dirname $file)"/"$Task"/"$OutPrefix"
  echo "$OutPrefix"
  mkdir -p $OutDir
  ExpectedOutput="$OutDir"/${OutPrefix}.ITS1.fasta

  Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  while [ "$Jobs" -gt 3 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -h -u theaven -n "$Task" | wc -l)
  done

  if [ ! -s "$ExpectedOutput" ]; then
    jobid=$(sbatch --job-name="$Task" --parsable ~/git_repos/Wrappers/unibz/run_ITSx.sh "$file" "$OutPrefix" "$OutDir")
    printf "%s\t%s\t "$Task" \t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> /home/clusterusers/theaven/slurm_log.tsv
  else
    echo "For $ID found: $ExpectedOutput" 
  fi
done
```

## Illumina sequencing

### Collecting data

Raw sequencing reads were retreived from the archive folder \\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\Illumina_16S_ITS and uploaded to the HPC:
```bash
ls /data/users/theaven/Ips_jam_project/raw_data/illumina/

for file in $(ls /data/users/theaven/Ips_jam_project/raw_data/illumina/*/*.fastq.gz); do
ID=$(basename $file | cut -d '_' -f1)
mkdir $(dirname $file)/$ID
mv $file $(dirname $file)/$ID/.
done
```
Sample metadata is here: "\\share.unibz.it\AppliedMolecularEntomologyLab\ips_typographus\Illumina_16S_ITS\Illumina_Ips_sequencing_list.xlsx"

## Quality Control and ASV Inference 

### FastQC 
The raw sequence reads were subjected to a quality control check using FastQC.
```bash
module load anaconda3
conda activate seqkit-2.10
seqkit stats /data/users/theaven/Ips_jam_project/raw_data/illumina/*/*/*fastq.gz

for ReadDir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/illumina/*/*); do
  Task=FastQC
  ID=$(echo "$ReadDir" | cut -d '/' -f9 )
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
### Cutadapt  
Primers were removed from the reads where present using Cutadapt. Primers used by Macrogen for the 16S V3-V4 and ITS3-ITS4 region are given at https://www.macrogen-europe.com/service/metagenome-sequencing

NOTE:The reads are a mix of paired and single end samples.
```bash
screen -r melanoneura
for ReadDir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/illumina/16S/*); do
  Task=CutAdapt
  ID=$(echo "$ReadDir" | cut -d '/' -f9 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
  OutDir="$(echo "$ReadDir" | sed 's@raw_data@qc_data@g')/"$Task""
  Forward_Primer=CCTACGGGNGGCWGCAG
  Reverse_Primer=GACTACHVGGGTATCTAATCC
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

for ReadDir in $(ls -d /data/users/theaven/Ips_jam_project/raw_data/illumina/ITS/*); do
  Task=CutAdapt
  ID=$(echo "$ReadDir" | cut -d '/' -f9 | sed 's@/@_@g')
    Reads=("$ReadDir"/*.fastq.gz)
  OutDir="$(echo "$ReadDir" | sed 's@raw_data@qc_data@g')/"$Task""
  Forward_Primer=GCATCGATGAAGAACGCAGC
  Reverse_Primer=TCCTCCGCTTATTGATATGC
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

seqkit stats /data/users/theaven/Ips_jam_project/qc_data/illumina/*/*/CutAdapt/*fastq.gz
```
### DADA2  <a name="7"></a>
The denoiser tool DADA2 was run to model and remove error patterns from the Illumina data: read ends where quality drops were trimmed, reads with high expected error (EE), or above a max EE threshold, or with any ambiguous bases (N) were discarded. DADA2’s removeBimeraDenovo function was also used to remove chimera from the ASV table. For paired reads DADA2 also merges forward and reverse reads in the overlapping region. When F and R disagree, the merge algorithm uses quality scores to pick the most likely base or discards the read. The output of DADA2 is an Alternative Sequence Variant (ASV) table as well as filtered denoised read fastq files.

DADA2 is an R package, the relavent files were downloaded:
```bash
for Dir in $(ls -d /data/users/theaven/Ips_jam_project/qc_data/illumina/*/*/CutAdapt); do
    if ls "$Dir"/*_1.trim.fastq.gz 1> /dev/null 2>&1 && ls "$Dir"/*_2.trim.fastq.gz 1> /dev/null 2>&1; then
        Out=/data/users/theaven/download_20260529/paired/$(echo "$Dir" | cut -d '/' -f9)
        mkdir -p "$Out"
        cp "$Dir"/*.fastq.gz "$Out"/.
    else
        Out=/data/users/theaven/download_20260319/single/$(echo "$Dir" | cut -d '/' -f9)
        mkdir -p "$Out"
        cp "$Dir"/*.fastq.gz "$Out"/.
    fi
done
```
"/data/users/theaven/download_20260529" was subseqeuntly downloaded to "C:\Users\THeaven\OneDrive - Scientific Network South Tyrol\R"

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

#Example reads plotted:
plotQualityProfile("download_20260529/paired/SHBB01881-1/SHBB01881-1_16S_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01881-1/SHBB01881-1_16S_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01887-2/SHBB01887-2_16S_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01887-2/SHBB01887-2_16S_2.trim.fastq.gz")  

plotQualityProfile("download_20260529/paired/SHBB01891-1/SHBB01891-1_16S_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01891-1/SHBB01891-1_16S_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01909-1/SHBB01909-1_16S_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01909-1/SHBB01909-1_16S_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01888-3/SHBB01888-3_16S_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01888-3/SHBB01888-3_16S_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01882-1/SHBB01882-1_16S_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01882-1/SHBB01882-1_16S_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB0-ve/SHBB0-ve_16S_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB0-ve/SHBB0-ve_16S_2.trim.fastq.gz") 


plotQualityProfile("download_20260529/paired/SHBB01881-1/SHBB01881-1_ITS_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01881-1/SHBB01881-1_ITS_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01887-2/SHBB01887-2_ITS_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01887-2/SHBB01887-2_ITS_2.trim.fastq.gz")  

plotQualityProfile("download_20260529/paired/SHBB01891-1/SHBB01891-1_ITS_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01891-1/SHBB01891-1_ITS_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01909-1/SHBB01909-1_ITS_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01909-1/SHBB01909-1_ITS_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01888-3/SHBB01888-3_ITS_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01888-3/SHBB01888-3_ITS_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB01882-1/SHBB01882-1_ITS_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB01882-1/SHBB01882-1_ITS_2.trim.fastq.gz") 

plotQualityProfile("download_20260529/paired/SHBB0-ve/SHBB0-ve_ITS_1.trim.fastq.gz") 
plotQualityProfile("download_20260529/paired/SHBB0-ve/SHBB0-ve_ITS_2.trim.fastq.gz") 
```

DADA2 quality profile plots summarise per-cycle base quality scores across reads. The solid orange line is median quality score at each base position, the solid turquoise line is mean quality score at each position, orange dashed lines show 10th and 90th percentiles. Q20 = ~1% error rate. Q30 = ~ 0.1% error rate.

Both mean (turquoise) and median (orange) are well above Q30 at all positions for every example plotted. The mean of reverse reads is slightly lower than forward reads, reverse reads are expected to be lower quality than forward reads.  

![Representative DADA2 quality profile plot SHBB01881-1 ITS reverse](figures/SHBB01881-1_ITS_2.png)

**Collect inputs to run DADA2 for 16S:**
```R
paired <- list.files(path = "download_20260529/paired", full.names = TRUE, recursive = TRUE)

paired_1 <- paired[grepl("\\_16S_1\\.trim.fastq.gz$", paired)]
paired_2 <- paired[grepl("\\_16S_2\\.trim.fastq.gz$", paired)]

get_samplename <- function(x) basename(dirname(x))

snF <- vapply(paired_1, get_samplename, character(1))
snR <- vapply(paired_2, get_samplename, character(1))
paired_samples <- intersect(snF, snR)
paired_lookup_table <- data.frame(
  sample = paired_samples,
  fnF = paired_1[match(paired_samples, snF)],
  fnR = paired_2[match(paired_samples, snR)],
  stringsAsFactors = FALSE
)
```
**Run DADA2 for 16S:**
As the reads appear to have very good quality strict settings can be used, truncLen was set to approximate read length post trimming.
```R
#Parameters
truncLen <- c(280, 280)  #Truncate reads to this length, then remove entirely
truncQ <- 20     #Truncate at first base with quality <= truncQ (higher = more conservative | lower = more permissive)
maxEE <- c(1,1)  #Maximum expected errors (4 for the reverse reads as these are known/expected to be lower quality)
maxN <- 0        #Maximum allowed N bases
rm.phix <- TRUE  #Remove PhiX reads (bacteriophage used as control in Illumina sequencing runs)
pool <- FALSE    #Pool samples for error rate estimation (FALSE = error model is learned per sample, this is more conservative and faster/less memory | pseudo = max sensitivity)
threads <- TRUE  #Use multithreading
outdir <- "download_20260529/results/dada2_output_16S"
in_df <- paired_lookup_table

#Create output directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Filter and trim reads (paired)
filtFs <- file.path(outdir, paste0(in_df$sample, "_F_filt.fastq.gz"))
filtRs <- file.path(outdir, paste0(in_df$sample, "_R_filt.fastq.gz"))
out2 <- filterAndTrim(fwd = in_df$fnF, filt = filtFs, rev = in_df$fnR, filt.rev = filtRs, truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, rm.phix = rm.phix, compress = TRUE, multithread = threads)

#Learn error rates, dereplicate, run DADA2 algorithm, and merge pairs
#Learn errors separately for F and R
errF <- learnErrors(filtFs, nbases = 5e10, multithread = threads)
errR <- learnErrors(filtRs, nbases = 5e10, multithread = threads)

#299554080 total bases in 1069836 reads from 25 samples will be used for learning the error rates. - 280,280,1,1 - 72%

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

#Plot merged read lengths after chimera removal:
asv_lengths <- nchar(colnames(seqtab.nochim))
summary(asv_lengths)
hist(asv_lengths, breaks=50, main="ASV lengths after chimera removal", xlab="Length (bp)")

# Save results
saveRDS(seqtab.nochim, file = file.path(outdir, "seqtab.nochim.rds"))
write.table(seqtab.nochim, file.path(outdir, "seqtab.nochim.tsv"), sep="\t", quote=FALSE, col.names=NA)

getN <- function(x) sum(getUniques(x))
track <- cbind(input    = out2[, "reads.in"], filtered = out2[, "reads.out"], denoisedF = sapply(dadaF, getN), denoisedR = sapply(dadaR, getN), merged = sapply(mergers, function(x) sum(x$abundance)))
rownames(track) <- in_df$sample
write.table(track, file.path(outdir, "read_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)

write.table(chim_track, file.path(outdir, "chimera_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
```
Remove ASV if is not at least 1% of reads in a sample or 0.1% of reads in multiple samples. (Sequencing errors and tag-jumps almost always show up in one sample only)
```R
ASV <- readRDS("download_20260529/results/dada2_output_16S/seqtab.nochim.rds")

prune_single_sample_low_abundance <- function(seqtab, min_rel = 0.01) {
  prevalence <- colSums(seqtab > 0) #count how many samples each ASV appears in
  sample_totals <- rowSums(seqtab) #find total reads per sample
  rel_abund <- sweep(seqtab, 1, sample_totals, "/") #ASV_reads / total_reads_in_sample
  rel_abund[is.na(rel_abund)] <- 0 #Set relative abundance to 0 where undefined
  max_rel <- apply(rel_abund, 2, max) #Finds the highest relative abundance ASV ever reaches (for single sample ASVs there is only 1 value)
  keep <- (max_rel >= 0.01) | (prevalence >= 2 & max_rel >= 0.001) #remove ASV if is not at least 1% of reads in a sample or 0.1% of reads in multiple samples
  seqtab[, keep, drop = FALSE]
}

ASV_pruned <- prune_single_sample_low_abundance(ASV, min_rel = 0.01)

ncol(ASV) #2291
ncol(ASV_pruned) #381
mean(rowSums(ASV_pruned > 0)) #63.32

dir.create("download_20260529/ASVs", showWarnings = FALSE, recursive = TRUE)

saveRDS(ASV_pruned, file = file.path("download_20260529/ASVs/ASV_asv_16s.rds"))
write.table(ASV_pruned, file.path("download_20260529/ASVs/ASV_asv_16s.tsv"), sep="\t", quote=FALSE, col.names=NA)

# seqtab_filt: samples x ASVs
asv_seqs <- colnames(ASV_pruned)
asv_headers <- paste0("ASV", seq_along(asv_seqs))
dna <- Biostrings::DNAStringSet(asv_seqs)
names(dna) <- asv_headers
Biostrings::writeXStringSet(dna, "download_20260529/ASVs/ASVs_16s.fasta")
write.csv(data.frame(ASV=asv_headers, Sequence=asv_seqs), "download_20260529/ASVs/ASV_16s_id_map.csv", row.names = FALSE)
```
**Collect inputs to run DADA2 for ITS:**
```R
paired <- list.files(path = "download_20260529/paired", full.names = TRUE, recursive = TRUE)

paired_1 <- paired[grepl("\\_ITS_1\\.trim.fastq.gz$", paired)]
paired_2 <- paired[grepl("\\_ITS_2\\.trim.fastq.gz$", paired)]

get_samplename <- function(x) basename(dirname(x))

snF <- vapply(paired_1, get_samplename, character(1))
snR <- vapply(paired_2, get_samplename, character(1))
paired_samples <- intersect(snF, snR)
paired_lookup_table <- data.frame(
  sample = paired_samples,
  fnF = paired_1[match(paired_samples, snF)],
  fnR = paired_2[match(paired_samples, snR)],
  stringsAsFactors = FALSE
)
```
**Run DADA2 for ITS:**
As the reads appear to have very good quality strict settings can be used, truncLen was set to approximate read length post trimming.
```R
#Parameters
truncLen <- c(280, 280)  #Truncate reads to this length, then remove entirely
truncQ <- 20     #Truncate at first base with quality <= truncQ (higher = more conservative | lower = more permissive)
maxEE <- c(1,1)  #Maximum expected errors (4 for the reverse reads as these are known/expected to be lower quality)
maxN <- 0        #Maximum allowed N bases
rm.phix <- TRUE  #Remove PhiX reads (bacteriophage used as control in Illumina sequencing runs)
pool <- FALSE    #Pool samples for error rate estimation (FALSE = error model is learned per sample, this is more conservative and faster/less memory | pseudo = max sensitivity)
threads <- TRUE  #Use multithreading
outdir <- "download_20260529/results/dada2_output_ITS"
in_df <- paired_lookup_table

#Create output directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Filter and trim reads (paired)
filtFs <- file.path(outdir, paste0(in_df$sample, "_F_filt.fastq.gz"))
filtRs <- file.path(outdir, paste0(in_df$sample, "_R_filt.fastq.gz"))
out2 <- filterAndTrim(fwd = in_df$fnF, filt = filtFs, rev = in_df$fnR, filt.rev = filtRs, truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, rm.phix = rm.phix, compress = TRUE, multithread = threads)

#Learn error rates, dereplicate, run DADA2 algorithm, and merge pairs
#Learn errors separately for F and R
errF <- learnErrors(filtFs, nbases = 5e10, multithread = threads)
errR <- learnErrors(filtRs, nbases = 5e10, multithread = threads)

#357098560 total bases in 1275352 reads from 25 samples will be used for learning the error rates. - 280,280,1,1 - 77%

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

#Plot merged read lengths after chimera removal:
asv_lengths <- nchar(colnames(seqtab.nochim))
summary(asv_lengths)
hist(asv_lengths, breaks=50, main="ASV lengths after chimera removal", xlab="Length (bp)")

# Save results
saveRDS(seqtab.nochim, file = file.path(outdir, "seqtab.nochim.rds"))
write.table(seqtab.nochim, file.path(outdir, "seqtab.nochim.tsv"), sep="\t", quote=FALSE, col.names=NA)

getN <- function(x) sum(getUniques(x))
track <- cbind(input    = out2[, "reads.in"], filtered = out2[, "reads.out"], denoisedF = sapply(dadaF, getN), denoisedR = sapply(dadaR, getN), merged = sapply(mergers, function(x) sum(x$abundance)))
rownames(track) <- in_df$sample
write.table(track, file.path(outdir, "read_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)

write.table(chim_track, file.path(outdir, "chimera_tracking.tsv"), sep="\t", quote=FALSE, col.names=NA)
```
Remove ASV if is not at least 1% of reads in a sample or 0.1% of reads in multiple samples. (Sequencing errors and tag-jumps almost always show up in one sample only)
```R
ASV <- readRDS("download_20260529/results/dada2_output_ITS/seqtab.nochim.rds")

prune_single_sample_low_abundance <- function(seqtab, min_rel = 0.01) {
  prevalence <- colSums(seqtab > 0) #count how many samples each ASV appears in
  sample_totals <- rowSums(seqtab) #find total reads per sample
  rel_abund <- sweep(seqtab, 1, sample_totals, "/") #ASV_reads / total_reads_in_sample
  rel_abund[is.na(rel_abund)] <- 0 #Set relative abundance to 0 where undefined
  max_rel <- apply(rel_abund, 2, max) #Finds the highest relative abundance ASV ever reaches (for single sample ASVs there is only 1 value)
  keep <- (max_rel >= 0.01) | (prevalence >= 2 & max_rel >= 0.001) #remove ASV if is not at least 1% of reads in a sample or 0.1% of reads in multiple samples
  seqtab[, keep, drop = FALSE]
}

ASV_pruned <- prune_single_sample_low_abundance(ASV, min_rel = 0.01)

ncol(ASV) #1099
ncol(ASV_pruned) #213
mean(rowSums(ASV_pruned > 0)) #50.52

dir.create("download_20260529/ASVs", showWarnings = FALSE, recursive = TRUE)

saveRDS(ASV_pruned, file = file.path("download_20260529/ASVs/ASV_asv_its.rds"))
write.table(ASV_pruned, file.path("download_20260529/ASVs/ASV_asv_its.tsv"), sep="\t", quote=FALSE, col.names=NA)

# seqtab_filt: samples x ASVs
asv_seqs <- colnames(ASV_pruned)
asv_headers <- paste0("ASV", seq_along(asv_seqs))
dna <- Biostrings::DNAStringSet(asv_seqs)
names(dna) <- asv_headers
Biostrings::writeXStringSet(dna, "download_20260529/ASVs/ASVs_its.fasta")
write.csv(data.frame(ASV=asv_headers, Sequence=asv_seqs), "download_20260529/ASVs/ASV_its_id_map.csv", row.names = FALSE)
```
### QIIME - plot
Plot rarefaction curve with qiime

Prepare inputs DADA2(R) -> QIIME2(HPC):

**16S**
```bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

ASV_dir=/data/users/theaven/Ips_jam_project/asvs/ASVs

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/make_qiime2_inputs.py \
  --counts "$ASV_dir"/ASV_asv_16s.tsv \
  --map "$ASV_dir"/ASV_16s_id_map.csv \
  --out-dir "$ASV_dir"/qiime_inputs \
  --out-prefix qiime_inputs_16s_

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i "$ASV_dir"/qiime_inputs/qiime_inputs_16s_ASV_table.tsv \
  -o "$ASV_dir"/qiime_inputs/ASV_table_16s.biom \
  --table-type="OTU table" \
  --to-hdf5

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path "$ASV_dir"/qiime_inputs/ASV_table_16s.biom \
  --output-path "$ASV_dir"/qiime_inputs/table_16s.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$ASV_dir"/ASVs_16s.fasta \
  --output-path "$ASV_dir"/rep-seqs_16s.qza

#Input metadata:
echo -e "#SampleID" > "$ASV_dir"/sample-metadata.tsv
head -n 1 "$ASV_dir"/qiime_inputs/qiime_inputs_16s_ASV_table.tsv | cut -f2- | tr '\t' '\n' >> "$ASV_dir"/sample-metadata.tsv
#Download, input metadata, upload

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_dir"/rep-seqs_16s.qza \
  --output-path "$ASV_dir"/rep-seqs_16s-fasta
```
Plot

```bash
#Build tree
srun -p bioagri  -c 4 --mem 64G --pty bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn
ASV_dir=/data/users/theaven/Ips_jam_project/asvs/ASVs
mkdir -p "$ASV_dir"/tmp
TMPDIR="$ASV_dir"/tmp \
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif \
qiime alignment mafft \
  --i-sequences "$ASV_dir"/rep-seqs_16s.qza \
  --o-alignment "$ASV_dir"/aligned-rep-seqs_16s.qza \
  --p-n-threads 4

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime alignment mask \
  --i-alignment "$ASV_dir"/aligned-rep-seqs_16s.qza \
  --o-masked-alignment "$ASV_dir"/masked-aligned-rep-seqs_16s.qza 

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime phylogeny fasttree \
  --i-alignment "$ASV_dir"/masked-aligned-rep-seqs_16s.qza \
  --o-tree "$ASV_dir"/unrooted-tree_16s.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime phylogeny midpoint-root \
  --i-tree "$ASV_dir"/unrooted-tree_16s.qza \
  --o-rooted-tree "$ASV_dir"/rooted-tree_16s.qza

#Find sequencing depth cuttoff
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table summarize \
  --i-table "$ASV_dir"/qiime_inputs/table_16s.qza \
  --o-visualization "$ASV_dir"/table_16s-summary.qzv \
  --m-sample-metadata-file "$ASV_dir"/sample-metadata.tsv
#0 samples are below 10,000 depth
#1 samples are below 20,000
#5 below 30,000
#14 below 40,000

#Plot rarefaction curve
awk 'NR>1 {for(i=2;i<=NF;i++) if($i>max) max=$i} END{print max}' "$ASV_dir"/ASV_asv_16s.tsv #43367
awk 'NR==1 {next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1, sum}' "$ASV_dir"/ASV_asv_16s.tsv | sort -k2 -n | sort -k2 -n | sort -k2 -n | head -3 #SHBB0-ve 15112, SHBB01898-1 22723, SHBB01915-1 24990
awk 'NR==1 {next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1, sum}' "$ASV_dir"/ASV_asv_16s.tsv | sort -k2 -n | sort -k2 -n | sort -k2 -n | tail -1 #SHBB01888-1 63031

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime diversity alpha-rarefaction \
  --i-table "$ASV_dir"/qiime_inputs/table_16s.qza \
  --i-phylogeny "$ASV_dir"/rooted-tree_16s.qza \
  --p-max-depth 15000 \
  --m-metadata-file "$ASV_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_dir"/alpha-rarefaction_16s.qzv 

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime diversity alpha-rarefaction \
  --i-table "$ASV_dir"/qiime_inputs/table_16s.qza \
  --i-phylogeny "$ASV_dir"/rooted-tree_16s.qza \
  --p-max-depth 30000 \
  --m-metadata-file "$ASV_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_dir"/alpha-rarefaction-30000_16s.qzv 
```
![Rarefaction curves for ITS seqeunces](figures/rarefaction_16s.png)

**ITS**
```bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

ASV_dir=/data/users/theaven/Ips_jam_project/asvs/ASVs

apptainer exec --bind /data:/data --bind /home/clusterusers/theaven:/home/clusterusers/theaven ~/git_repos/Containers/python3.sif python ~/git_repos/Scripts/unibz/make_qiime2_inputs.py \
  --counts "$ASV_dir"/ASV_asv_its.tsv \
  --map "$ASV_dir"/ASV_its_id_map.csv \
  --out-dir "$ASV_dir"/qiime_inputs \
  --out-prefix qiime_inputs_its_

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i "$ASV_dir"/qiime_inputs/qiime_inputs_its_ASV_table.tsv \
  -o "$ASV_dir"/qiime_inputs/ASV_table_its.biom \
  --table-type="OTU table" \
  --to-hdf5

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path "$ASV_dir"/qiime_inputs/ASV_table_its.biom \
  --output-path "$ASV_dir"/qiime_inputs/table_its.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$ASV_dir"/ASVs_its.fasta \
  --output-path "$ASV_dir"/rep-seqs_its.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_dir"/rep-seqs_its.qza \
  --output-path "$ASV_dir"/rep-seqs_its-fasta
```
Plot

```bash
#Build tree
srun -p bioagri  -c 4 --mem 64G --pty bash
module load apptainer/1.4.1-gcc-13.3.0-3coysxn
ASV_dir=/data/users/theaven/Ips_jam_project/asvs/ASVs
mkdir -p "$ASV_dir"/tmp
TMPDIR="$ASV_dir"/tmp \
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif \
qiime alignment mafft \
  --i-sequences "$ASV_dir"/rep-seqs_its.qza \
  --o-alignment "$ASV_dir"/aligned-rep-seqs_its.qza \
  --p-n-threads 4

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime alignment mask \
  --i-alignment "$ASV_dir"/aligned-rep-seqs_its.qza \
  --o-masked-alignment "$ASV_dir"/masked-aligned-rep-seqs_its.qza 

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime phylogeny fasttree \
  --i-alignment "$ASV_dir"/masked-aligned-rep-seqs_its.qza \
  --o-tree "$ASV_dir"/unrooted-tree_its.qza

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime phylogeny midpoint-root \
  --i-tree "$ASV_dir"/unrooted-tree_its.qza \
  --o-rooted-tree "$ASV_dir"/rooted-tree_its.qza

#Find sequencing depth cuttoff
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime feature-table summarize \
  --i-table "$ASV_dir"/qiime_inputs/table_its.qza \
  --o-visualization "$ASV_dir"/table_its-summary.qzv \
  --m-sample-metadata-file "$ASV_dir"/sample-metadata.tsv
#0 samples are below 10,000 depth
#1 samples are below 20,000
#4 below 30,000
#8 below 40,000

#Plot rarefaction curve
awk 'NR>1 {for(i=2;i<=NF;i++) if($i>max) max=$i} END{print max}' "$ASV_dir"/ASV_asv_its.tsv #41027
awk 'NR==1 {next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1, sum}' "$ASV_dir"/ASV_asv_its.tsv | sort -k2 -n | sort -k2 -n | sort -k2 -n | head -3 #SHBB01895-1 16411, SHBB0-ve 20729, SHBB01915-1 21150
awk 'NR==1 {next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1, sum}' "$ASV_dir"/ASV_asv_its.tsv | sort -k2 -n | sort -k2 -n | sort -k2 -n | tail -1 #SHBB01898-1 63188

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime diversity alpha-rarefaction \
  --i-table "$ASV_dir"/qiime_inputs/table_its.qza \
  --i-phylogeny "$ASV_dir"/rooted-tree_its.qza \
  --p-max-depth 16400 \
  --m-metadata-file "$ASV_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_dir"/alpha-rarefaction_its.qzv 

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime diversity alpha-rarefaction \
  --i-table "$ASV_dir"/qiime_inputs/table_its.qza \
  --i-phylogeny "$ASV_dir"/rooted-tree_its.qza \
  --p-max-depth 30000 \
  --m-metadata-file "$ASV_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_dir"/alpha-rarefaction-30000_its.qzv 
```
![Rarefaction curves for ITS seqeunces](figures/rarefaction-its.png)

## Taxonomy Assignment <a name="9"></a>

Ready-made IDTAXA training set is publically available for SILVA SSU and was downloaded from https://www2.decipher.codes/Downloads.html - Accessed 07/06/2026 - "   SILVA SSU r138.2 (modified) " and "   UNITE 2025 (unmodified) "

### IDTAXA <a name="10"></a>

IDTAXA (DECIPHER package in R) - probabilistic sequence classification, model-based learning - usually more accurate and conservative than BLAST or naive Bayes. IDTAXA does not require trimming and handles full-length sequences correctly. IDTAXA is slower than QIIME2 NB, but usually much more accurate. If trained on full-length, IDTAXA often backs off to a safer rank rather than confidently guessing too deep.

NOTE: the 'species' data in SILVA can be questionable as users list host organism in this field or write things like 'metagenome'. We will need to be carefull when interpreting the Assignment results at the species level. Also, for some entries species is given but not all the higher level taxa, this can cause problems because IDTAXA likes unique species to have only one taxonomic path. Other entries like Incertae Sedis also had to be made unambiguous. I have therefore changed every column to give the full ':' seperated taxonomic path to that level in order to bypass this problem, this tranformation will presumably need to be reversed in the outputs of IDTAXA in order to make them readable. 

**16S**

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DECIPHER","stringr", "Biostrings"))

library(DECIPHER)
library(Biostrings)
library(stringr)

setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")

load("SILVA_SSU_r138_2_2024.RData")

SSU_trainingSet <- trainingSet

SSU_ASVs <- readDNAStringSet("download_20260529/ASVs/ASVs_16s.fasta")

SSU_tax_idtaxa <- IdTaxa(
  SSU_ASVs,
  SSU_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  #% of bootstraps supporting assignment, raise to 70–80 for fewer false positives
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
write.table(SSU_tax_tab, file = "download_20260529/ASVs/SSU_tax_tab_corrected_16s.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")

#convert to qiime format for plotting
tax_df <- data.frame(
  Feature.ID = names(SSU_tax_idtaxa),
  Taxon = sapply(SSU_tax_idtaxa, function(x) {
    # Skip the first rank (root)
    ranks <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    paste0(ranks, x$taxon[-1])[1:length(ranks)] |> paste(collapse = "; ")
  }),
  Confidence = sapply(SSU_tax_idtaxa, function(x) {
    min(x$confidence, na.rm = TRUE) / 100
  })
)

write.table(
  tax_df,
  "download_20260529/ASVs/idtaxa_taxonomy_16s.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```
Curate the taxonomy file: replace g_endosymbionts with g_unclassified_f, replace Incertae Sedis entries with x_unclassified_x, propogate x_unclassified_x entries down ranks. ->  idtaxa_taxonomy2.tsv

**ITS**

```R
load("UNITE_v2025.RData")

ITS_trainingSet <- trainingSet

ITS_ASVs <- readDNAStringSet("download_20260529/ASVs/ASVs_its.fasta")

ITS_tax_idtaxa <- IdTaxa(
  ITS_ASVs,
  ITS_trainingSet,
  strand = "both",
  bootstraps = 100 , #default - maximum number of bootstrap replicates to perform for each sequence
  processors = NULL, #Use all available cores
  threshold = 60 ,  #% of bootstraps supporting assignment, raise to 70–80 for fewer false positives
  verbose = TRUE 
)

ITS_tax_tab <- t(vapply(ITS_tax_idtaxa, extract_tax,
                    FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
ITS_tax_tab <- as.data.frame(ITS_tax_tab, stringsAsFactors = FALSE)
rownames(ITS_tax_tab) <- names(ITS_ASVs)
write.table(ITS_tax_tab, file = "download_20260529/ASVs/ITS_tax_tab_corrected_its.tsv", sep = "\t", row.names = TRUE, quote = FALSE, na = "NA")

#convert to qiime format for plotting
tax_df2 <- data.frame(
  Feature.ID = names(ITS_tax_idtaxa),
  Taxon = sapply(ITS_tax_idtaxa, function(x) {
    # Skip the first rank (root)
    ranks <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    paste0(ranks, x$taxon[-1])[1:length(ranks)] |> paste(collapse = "; ")
  }),
  Confidence = sapply(ITS_tax_idtaxa, function(x) {
    min(x$confidence, na.rm = TRUE) / 100
  })
)

write.table(
  tax_df2,
  "download_20260529/ASVs/idtaxa_taxonomy_its.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```
### Curate

Curate the taxonomy file: replace g_endosymbionts with g_unclassified_f, replace Incertae Sedis entries with x_unclassified_x, propogate x_unclassified_x entries down ranks. 
Data entry error for ASV370 - repeat
```bash
ASV_dir=/data/users/theaven/Ips_jam_project/asvs/ASVs
module load apptainer/1.4.1-gcc-13.3.0-3coysxn

#import
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path "$ASV_dir"/idtaxa_taxonomy_16s.tsv \
  --output-path "$ASV_dir"/idtaxa_taxonomy_16s.qza \
  --input-format HeaderlessTSVTaxonomyFormat

apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path "$ASV_dir"/idtaxa_taxonomy_its.tsv \
  --output-path "$ASV_dir"/idtaxa_taxonomy_its.qza \
  --input-format HeaderlessTSVTaxonomyFormat

#remove organelle and plant hits
apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa filter-table \
  --i-table "$ASV_dir"/qiime_inputs/table_16s.qza \
  --i-taxonomy "$ASV_dir"/idtaxa_taxonomy_16s.qza \
  --p-exclude mitochondria,chloroplast,Viridiplantae \
  --o-filtered-table "$ASV_dir"/table-no-mitochondria-chloroplast_16s.qza

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa filter-table \
  --i-table "$ASV_dir"/qiime_inputs/table_its.qza \
  --i-taxonomy "$ASV_dir"/idtaxa_taxonomy_its.qza \
  --p-exclude mitochondria,chloroplast,Viridiplantae \
  --o-filtered-table "$ASV_dir"/table-no-mitochondria-chloroplast_its.qza

#collapse to genus/family level:
apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa collapse \
  --i-table "$ASV_dir"/table-no-mitochondria-chloroplast_its.qza \
  --i-taxonomy "$ASV_dir"/idtaxa_taxonomy_its.qza \
  --p-level 6 \
  --o-collapsed-table genus-table_its.qza

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path genus-table_its.qza \
  --output-path genus_export_its

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i genus_export_its/feature-table.biom \
  -o genus-table_its.tsv \
  --to-tsv

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa collapse \
  --i-table "$ASV_dir"/table-no-mitochondria-chloroplast_its.qza \
  --i-taxonomy "$ASV_dir"/idtaxa_taxonomy_its.qza \
  --p-level 5 \
  --o-collapsed-table family-table_its.qza

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path family-table_its.qza \
  --output-path family_export_its

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i family_export_its/feature-table.biom \
  -o family-table_its.tsv \
  --to-tsv

```

### Plot

**16S**

```bash
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa barplot \
  --i-table "$ASV_dir"/table-no-mitochondria-chloroplast_16s.qza \
  --i-taxonomy "$ASV_dir"/idtaxa_taxonomy_16s.qza \
  --m-metadata-file "$ASV_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_dir"/idtaxa-barplot_16s.qzv
```
![Relative abundance plots for ITS seqeunces](figures/bar-16s.png)

Download to plot with R:
```bash
down_dir=/data/users/theaven/Ips_jam_project/down_20260609
mkdir "$down_dir"

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_dir"/table-no-mitochondria-chloroplast_16s.qza \
  --output-path "$ASV_dir"/exported-feature-table_16s

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i "$ASV_dir"/exported-feature-table_16s/feature-table.biom \
  -o "$down_dir"/feature-table_16s.tsv \
  --to-tsv

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_dir"/idtaxa_taxonomy_16s.qza \
  --output-path "$down_dir"/exported-taxonomy_16s

cp /data/users/theaven/Ips_jam_project/asvs/ASVs/sample-metadata.tsv "$down_dir"/

###

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_dir"/qiime_inputs/table_16s.qza \
  --output-path "$ASV_dir"/exported-feature-table_16sc

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i "$ASV_dir"/exported-feature-table_16sc/feature-table.biom \
  -o "$down_dir"/feature-table_16sc.tsv \
  --to-tsv
```

Abundance plot:
```R
install.packages(c("tidyverse"))
install.packages("BiocManager")
install.packages("ggrepel")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
install.packages("ComplexUpset")
install.packages("pheatmap")
library(ComplexUpset)
library(phyloseq)
library(tidyverse)
library(readr)
library(tibble)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)

setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")
set.seed(1)

# Load table
otu <- read.table("down_20260609/feature-table_16sc.tsv", header=TRUE, row.names=1, sep="\t", comment.char="")
otu <- as.matrix(otu)
colnames(otu) <- gsub("\\.", "-", colnames(otu))

# Load metadata
meta <- read_tsv(
  "down_20260609/sample-metadata.tsv",
  comment = "",  
  show_col_types = FALSE
)
meta <- column_to_rownames(meta, var = "#SampleID")

tax1 <- read.table("down_20260609/exported-taxonomy_16s/taxonomy.tsv", 
                  header = TRUE, 
                  sep = "\t", 
                  row.names = 1)
tax_split <- tax1 %>%
  separate(Taxon, 
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
           sep = ";", 
           fill = "right")

# Create objects
OTU <- otu_table(otu, taxa_are_rows=TRUE)
SAM <- sample_data(meta)
TAX <- tax_table(as.matrix(tax_split))
ps <- phyloseq(OTU, SAM, TAX)

tax_table(ps) <- apply(tax_table(ps), 2, trimws)
tax_table(ps)[, "Genus"] <- gsub("^s__", "g__", tax_table(ps)[, "Genus"])

ps_genus <- tax_glom(ps, taxrank = "Genus")

ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

df <- psmelt(ps_rel)

taxa_abund <- tapply(df$Abundance, df$Genus, sum)

top10 <- names(sort(taxa_abund, decreasing = TRUE))[1:10]

keep_taxa <- unique(c(top10))

df$Genus <- as.character(df$Genus)

df$Genus[!df$Genus %in% keep_taxa] <- "Other"

all_taxa <- unique(df$Genus)

auto_taxa <- setdiff(all_taxa, "Other")

auto_cols <- setNames(c(
  "darkblue",  "#008000", "sienna3", "orange",
  "deeppink4",  "skyblue3", "red3"
  , "maroon3", "wheat3", "yellow"
), auto_taxa)

other_col <- c("Other" = "grey80")

final_cols <- c(auto_cols, other_col)

df$Genus <- factor(df$Genus, levels = names(final_cols))

#remove sequencing blank from the plot
df_sub <- subset(df, treatment != "blank")

# create hierarchical grouping key
df_sub$treatment <- as.character(df_sub$treatment)
df_sub$Group <- paste(df_sub$treatment, df_sub$Sample, sep = " ")
df_sub$Label <- paste(df_sub$treatment, df_sub$Sample)

df_sub$Group <- factor(df_sub$Group, levels = unique(df_sub$Group))
df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

df_sub <- df_sub[order(df_sub$treatment, df_sub$Sample), ]
df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

ggplot(df_sub, aes(x = Label, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = final_cols) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
```
Alpha diversity:
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

meta_f <- meta[meta$treatment != "blank", , drop = FALSE]  
ps_f <- prune_samples(rownames(meta_f), ps)
identical(sample_names(ps_f), rownames(meta_f))
alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = treatment, y = Shannon, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = treatment),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  scale_fill_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")

#######

by(alpha_df$Shannon, alpha_df$treatment, shapiro.test)
#Control - p-value = 0.07682 - borderline
#Insecticide - p-value = 0.7613 - normal
#Microsap - p-value = 0.03346 - not normal

pairwise.wilcox.test(alpha_df$Shannon,
                     alpha_df$treatment,
                     p.adjust.method = "BH")
#No significant differences in Shannon diversity between any pair of treatments
#            Control Insecticide
#Insecticide 0.075   -          
#Microsap    0.382   0.075   
```
beta diversity:

```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    max.overlaps = Inf,
    force = 2
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

#######

betadisper_res <- betadisper(bray, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.00600 0.0029997 0.1817    999  0.818
#Residuals 21 0.34679 0.0165138
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(bray ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.6544 0.07565 0.8593  0.675
#Residual 21   7.9968 0.92435              
#Total    23   8.6512 1.00000    
#No significant difference in community composition between treatments, treatment explains only ~7.6% of variation in community composition     

#######

ps_simple <- phyloseq::phyloseq(
  phyloseq::otu_table(ps_f),
  phyloseq::sample_data(ps_f)
)
ps_pa <- transform_sample_counts(ps_simple, function(x) as.numeric(x > 0))
jaccard <- phyloseq::distance(ps_pa, method = "jaccard")
ord <- ordinate(ps_f, method = "PCoA", distance = jaccard)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")


  ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    max.overlaps = Inf,
    force = 2
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))
#######


betadisper_res <- betadisper(jaccard, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.005156 0.0025779 0.7559    999  0.476
#Residuals 21 0.071622 0.0034106
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(jaccard ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs    R2      F Pr(>F)  
#Model     2   0.9574 0.116 1.3779  0.013 *
#Residual 21   7.2961 0.884                
#Total    23   8.2535 1.000  
#Significant difference in community composition between treatments, treatment explains ~11.6% of variation     
```
Treatments affect which taxa are present/absent, not their relative abundances - community membership shifts but dominant taxa abundances are stable. Treatments do not strongly reshape dominant community structure, but they do cause taxa turnover (gains/losses of rarer taxa) - hidden signal in rare taxa.

Presence/abscence - bubble plot:
```R
#build matrix for Genus level
ps_f_genus <- tax_glom(ps_f, taxrank = "Genus")
ps_f_genus <- filter_taxa(ps_f_genus, function(x) sum(x) > 0, TRUE)

mat <- as(otu_table(ps_f_genus), "matrix")
if (taxa_are_rows(ps_f_genus)) {
  mat <- t(mat)
}

pa_mat <- (mat > 0) * 1

meta_mat <- meta_f[match(rownames(pa_mat), rownames(meta_f)), , drop = FALSE]

group <- meta_mat$treatment
groups <- unique(group)
stopifnot(all(rownames(pa_mat) == rownames(meta_mat)))

#Pairwise Fisher tests (ALL pairs)
pair_list <- combn(groups, 2, simplify = FALSE)

pairwise_p <- lapply(pair_list, function(grp) {
  g1 <- grp[1]
  g2 <- grp[2]
  idx <- group %in% c(g1, g2)
  apply(pa_mat[idx, , drop = FALSE], 2, function(x) {
    tab <- table(x, group[idx])
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA)
    if (all(tab == 0)) return(NA)
    fisher.test(tab)$p.value
  })
})

names(pairwise_p) <- sapply(pair_list, paste, collapse = "_vs_")

pairwise_p_adj <- lapply(pairwise_p, function(p) {
  p.adjust(p, method = "BH")
})

#Extract significant taxa (ANY comparison)
sig_taxa <- unique(unlist(lapply(pairwise_p_adj, function(p) {
  names(p)[which(p < 0.05 & !is.na(p))]
})))

prev_list <- lapply(groups, function(g) {
  colMeans(pa_mat[group == g, sig_taxa, drop = FALSE])
})
names(prev_list) <- groups

tax <- as.data.frame(tax_table(ps_f_genus))
tax_labels <- tax$Genus
names(tax_labels) <- rownames(tax)

df <- data.frame(
  Taxon = sig_taxa,
  Label = tax_labels[sig_taxa]
)

for (g in groups) {
  df[[paste0(g, "_prev")]] <- prev_list[[g]]
}

#Long format for plotting
df_long <- pivot_longer(
  df,
  cols = ends_with("_prev"),
  names_to = "Group",
  values_to = "Prevalence"
)

df_long$Group <- gsub("_prev", "", df_long$Group)

df_long$Label[is.na(df_long$Label)] <- df_long$Taxon

ggplot(df_long, aes(
  x = Group,
  y = Label,
  size = Prevalence,
  color = Group
)) +
  geom_point(alpha = 0.85) +
  scale_size(range = c(2, 10)) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Genus",
    size = "Prevalence",
    color = "Group",
    title = "Significant taxa (pairwise Fisher tests, FDR < 0.05)"
  )
```
No significant taxa - hard with only 8 reps

Presence/abscence + adundance - heatmap:
```R
tax_tab <- as.data.frame(tax_table(ps_f_genus))

tax_names <- tax_tab$Genus
tax_names[is.na(tax_names) | tax_names == ""] <- "Unknown"
tax_names <- make.unique(tax_names)

mat_asv <- mat
colnames(mat_asv) <- tax_names

mat_genus_rel <- sweep(mat_asv, 1, rowSums(mat_asv), "/")
mat_genus_rel_log <- log10(mat_genus_rel + 1e-6)

#Order samples
meta_mat <- meta[rownames(mat_genus_rel_log), , drop = FALSE ]

ord <- order(meta_mat$treatment)

mat_ordered <- mat_genus_rel_log[ord, ]
meta_ordered <- meta_mat[ord, , drop = FALSE  ]

# --- Define colors: 0 = white, then blue → red ---
# Avoid including 0 in gradient
nonzero_vals <- mat_ordered[mat_ordered > 0]

# --- Row annotations ---
annotation_row <- data.frame(
  Treatment = meta_ordered$treatment
)

rownames(annotation_row) <- rownames(mat_ordered)

#Optional: gaps
gaps <- cumsum(table(meta_ordered$treatment))

pheatmap(
  mat_ordered,
  color = colorRampPalette(c("white", "blue", "red"))(100),
  breaks = seq(
    min(mat_ordered, na.rm = TRUE),
    max(mat_ordered, na.rm = TRUE),
    length.out = 101
  ),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  gaps_row = gaps,
  border_color = "grey90",
  fontsize_col = 6
)
```

**ITS**

```bash
apptainer exec ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime taxa barplot \
  --i-table "$ASV_dir"/table-no-mitochondria-chloroplast_its.qza \
  --i-taxonomy "$ASV_dir"/idtaxa_taxonomy_its.qza \
  --m-metadata-file "$ASV_dir"/sample-metadata.tsv \
  --o-visualization "$ASV_dir"/idtaxa-barplot_its.qzv
```
![Relative abundance plots for ITS seqeunces](figures/bar-its.png)

Download to plot with R:
```bash
apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_dir"/table-no-mitochondria-chloroplast_its.qza \
  --output-path "$ASV_dir"/exported-feature-table_its

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif biom convert \
  -i "$ASV_dir"/exported-feature-table_its/feature-table.biom \
  -o "$down_dir"/feature-table_its.tsv \
  --to-tsv

apptainer exec --bind /data ~/git_repos/Containers/qiime2-amplicon-2025.7.sif qiime tools export \
  --input-path "$ASV_dir"/idtaxa_taxonomy_its.qza \
  --output-path "$down_dir"/exported-taxonomy_its
```
```R
setwd("C:/Users/THeaven/OneDrive - Scientific Network South Tyrol/R")
set.seed(1)

# Load table
otu <- read.table("down_20260609/feature-table_its.tsv", header=TRUE, row.names=1, sep="\t", comment.char="")
otu <- as.matrix(otu)
colnames(otu) <- gsub("\\.", "-", colnames(otu))

# Load metadata
meta <- read_tsv(
  "down_20260609/sample-metadata.tsv",
  comment = "",  
  show_col_types = FALSE
)
meta <- column_to_rownames(meta, var = "#SampleID")

tax1 <- read.table("down_20260609/exported-taxonomy_its/taxonomy.tsv", 
                  header = TRUE, 
                  sep = "\t", 
                  row.names = 1)
tax_split <- tax1 %>%
  separate(Taxon, 
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
           sep = ";", 
           fill = "right")

# Create objects
OTU <- otu_table(otu, taxa_are_rows=TRUE)
SAM <- sample_data(meta)
TAX <- tax_table(as.matrix(tax_split))
ps <- phyloseq(OTU, SAM, TAX)

tax_table(ps) <- apply(tax_table(ps), 2, trimws)
tax_table(ps)[, "Genus"] <- gsub("^s__", "g__", tax_table(ps)[, "Genus"])

ps_genus <- tax_glom(ps, taxrank = "Genus")

ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

df <- psmelt(ps_rel)

taxa_abund <- tapply(df$Abundance, df$Genus, sum)

top10 <- names(sort(taxa_abund, decreasing = TRUE))[1:10]

keep_taxa <- unique(c(top10))

df$Genus <- as.character(df$Genus)

df$Genus[!df$Genus %in% keep_taxa] <- "Other"

all_taxa <- unique(df$Genus)

auto_taxa <- setdiff(all_taxa, "Other")

auto_cols <- setNames(c(
  "darkblue",  "#008000", "sienna3", "orange",
  "deeppink4",  "skyblue3", "red3"
  , "maroon3", "wheat3", "yellow"
), auto_taxa)

other_col <- c("Other" = "grey80")

final_cols <- c(auto_cols, other_col)

df$Genus <- factor(df$Genus, levels = names(final_cols))

#remove sequencing blank from the plot
df_sub <- subset(df, treatment != "blank")

# create hierarchical grouping key
df_sub$treatment <- as.character(df_sub$treatment)
df_sub$Group <- paste(df_sub$treatment, df_sub$Sample, sep = " ")
df_sub$Label <- paste(df_sub$treatment, df_sub$Sample)

df_sub$Group <- factor(df_sub$Group, levels = unique(df_sub$Group))
df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

df_sub <- df_sub[order(df_sub$treatment, df_sub$Sample), ]
df_sub$Label <- factor(df_sub$Label, levels = unique(df_sub$Label))

ggplot(df_sub, aes(x = Label, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = final_cols) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
```
Alpha diversity:
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

meta_f <- meta[meta$treatment != "blank", , drop = FALSE]  
ps_f <- prune_samples(rownames(meta_f), ps)
identical(sample_names(ps_f), rownames(meta_f))
alpha_df <- get_alpha(ps_f, meta_f)

ggplot(alpha_df, aes(x = treatment, y = Shannon, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(aes(color = treatment),
              width = 0.15,
              alpha = 0.7,
              size = 2) +
  scale_fill_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  guides(color = "none")

#######

by(alpha_df$Shannon, alpha_df$treatment, shapiro.test)
#Control - p-value = 0.7432 - normal
#Insecticide - p-value = 0.7337 - normal
#Microsap - p-value = 0.09381 - normal

pairwise.t.test(alpha_df$Shannon,
                alpha_df$treatment,
                p.adjust.method = "BH",
                pool.sd = FALSE)
#            Control Insecticide
#Insecticide 0.88    -          
#Microsap    0.88    0.92

pairwise.wilcox.test(alpha_df$Shannon,
                     alpha_df$treatment,
                     p.adjust.method = "BH")
#No significant differences in Shannon diversity between any pair of treatments
#            Control Insecticide
#Insecticide 0.96    -          
#Microsap    0.96    0.96 
```
beta diversity:

```R
bray <- phyloseq::distance(ps_f, method = "bray")
ord <- ordinate(ps_f, method = "PCoA", distance = bray)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    max.overlaps = Inf,
    force = 2
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5))

#######

betadisper_res <- betadisper(bray, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
#Groups     2 0.04012 0.020058 1.019    999  0.357
#Residuals 21 0.41336 0.019684
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(bray ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.8518 0.10682 1.2558  0.131
#Residual 21   7.1218 0.89318              
#Total    23   7.9736 1.00000    
#No significant difference in community composition between treatments, treatment explains only ~10.7% of variation in community composition     

#######

ps_simple <- phyloseq::phyloseq(
  phyloseq::otu_table(ps_f),
  phyloseq::sample_data(ps_f)
)
ps_pa <- transform_sample_counts(ps_simple, function(x) as.numeric(x > 0))
jaccard <- phyloseq::distance(ps_pa, method = "jaccard")
ord <- ordinate(ps_f, method = "PCoA", distance = jaccard)
pcoa_df <- as.data.frame(ord$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$Sample <- rownames(pcoa_df)
meta_f$Sample <- rownames(meta_f)
pcoa_df <- merge(pcoa_df, meta_f, by = "Sample")


  ggplot(pcoa_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Control" = "royalblue",
    "Insecticide" = "#CC6666",
    "Microsap" = "#66CC66"
  )) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )

#######


betadisper_res <- betadisper(jaccard, meta_f$treatment)
permutest(betadisper_res, permutations = 999)
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.017639 0.0088195 1.3242    999  0.277
#Residuals 21 0.139863 0.0066601
#No evidence of differences in dispersion - no difference in within-group variability
boxplot(betadisper_res)
adonis2(jaccard ~ treatment, data = meta_f, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.5967 0.10077 1.1766  0.163
#Residual 21   5.3251 0.89923              
#Total    23   5.9218 1.00000  
#No significant difference in community composition between treatments, treatment explains ~10.1% of variation     
```
No differences

Presence/abscence - bubble plot:
```R
#build matrix for Genus level
ps_f_genus <- tax_glom(ps_f, taxrank = "Genus")
ps_f_genus <- filter_taxa(ps_f_genus, function(x) sum(x) > 0, TRUE)

mat <- as(otu_table(ps_f_genus), "matrix")
if (taxa_are_rows(ps_f_genus)) {
  mat <- t(mat)
}

pa_mat <- (mat > 0) * 1

meta_mat <- meta_f[match(rownames(pa_mat), rownames(meta_f)), , drop = FALSE]

group <- meta_mat$treatment
groups <- unique(group)
stopifnot(all(rownames(pa_mat) == rownames(meta_mat)))

#Pairwise Fisher tests (ALL pairs)
pair_list <- combn(groups, 2, simplify = FALSE)

pairwise_p <- lapply(pair_list, function(grp) {
  g1 <- grp[1]
  g2 <- grp[2]
  idx <- group %in% c(g1, g2)
  apply(pa_mat[idx, , drop = FALSE], 2, function(x) {
    tab <- table(x, group[idx])
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA)
    if (all(tab == 0)) return(NA)
    fisher.test(tab)$p.value
  })
})

names(pairwise_p) <- sapply(pair_list, paste, collapse = "_vs_")

pairwise_p_adj <- lapply(pairwise_p, function(p) {
  p.adjust(p, method = "BH")
})

#Extract significant taxa (ANY comparison)
sig_taxa <- unique(unlist(lapply(pairwise_p_adj, function(p) {
  names(p)[which(p < 0.05 & !is.na(p))]
})))

prev_list <- lapply(groups, function(g) {
  colMeans(pa_mat[group == g, sig_taxa, drop = FALSE])
})
names(prev_list) <- groups

tax <- as.data.frame(tax_table(ps_f_genus))
tax_labels <- tax$Genus
names(tax_labels) <- rownames(tax)

df <- data.frame(
  Taxon = sig_taxa,
  Label = tax_labels[sig_taxa]
)

for (g in groups) {
  df[[paste0(g, "_prev")]] <- prev_list[[g]]
}

#Long format for plotting
df_long <- pivot_longer(
  df,
  cols = ends_with("_prev"),
  names_to = "Group",
  values_to = "Prevalence"
)

df_long$Group <- gsub("_prev", "", df_long$Group)

df_long$Label[is.na(df_long$Label)] <- df_long$Taxon

ggplot(df_long, aes(
  x = Group,
  y = Label,
  size = Prevalence,
  color = Group
)) +
  geom_point(alpha = 0.85) +
  scale_size(range = c(2, 10)) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Genus",
    size = "Prevalence",
    color = "Group",
    title = "Significant taxa (pairwise Fisher tests, FDR < 0.05)"
  )
```
No significant taxa - hard with only 8 reps

Presence/abscence + adundance - heatmap:
```R
tax_tab <- as.data.frame(tax_table(ps_f_genus))

tax_names <- tax_tab$Genus
tax_names[is.na(tax_names) | tax_names == ""] <- "Unknown"
tax_names <- make.unique(tax_names)

mat_asv <- mat
colnames(mat_asv) <- tax_names

mat_genus_rel <- sweep(mat_asv, 1, rowSums(mat_asv), "/")
mat_genus_rel_log <- log10(mat_genus_rel + 1e-6)

#Order samples
meta_mat <- meta[rownames(mat_genus_rel_log), , drop = FALSE ]

ord <- order(meta_mat$treatment)

mat_ordered <- mat_genus_rel_log[ord, ]
meta_ordered <- meta_mat[ord, , drop = FALSE  ]

# --- Define colors: 0 = white, then blue → red ---
# Avoid including 0 in gradient
nonzero_vals <- mat_ordered[mat_ordered > 0]

# --- Row annotations ---
annotation_row <- data.frame(
  Treatment = meta_ordered$treatment
)

rownames(annotation_row) <- rownames(mat_ordered)

#Optional: gaps
gaps <- cumsum(table(meta_ordered$treatment))

pheatmap(
  mat_ordered,
  color = colorRampPalette(c("white", "blue", "red"))(100),
  breaks = seq(
    min(mat_ordered, na.rm = TRUE),
    max(mat_ordered, na.rm = TRUE),
    length.out = 101
  ),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  gaps_row = gaps,
  border_color = "grey90",
  fontsize_col = 7
)
```