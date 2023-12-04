## PKINGS Terra Processing
## July 2022

## Copy Data to Google Cloud

All files stored on Amazon S3 need to be transferred to Google Cloud.
1. Setup Terra.bio Workspace and obtain bucket for work
2. Run the following:
  inputdata='s3_file_payload.txt'
  bucket='fc-secure-001a4fd5-e796-46ba-8aa8-de834384e0bd'
  sbatch -a 0-500 --export=inputdata=$inputdata,bucket=$bucket submit_copy.sb


## Configure Reference Data:

### Create Reference

For this run, I ordered the pkings_0.2 assembly using the Brienomyrus brachyistius chromsome-level assembly.  
For details see github repo:
`ragtag_2022`

### Create Interval files

#### setup virtual machines
docker run -d -it --name gatk --mount type=bind,source=/Users/jasongallant/Desktop/for_github/ragtag_2022/pk_0.2_to_bbrach_final/,target=/jrgdata broadinstitute/gatk:4.0.1.0

#### prepare files

docker exec gatk cat /jrgdata/ragtag.scaffold.fasta.fai | cut -f1 > /Users/jasongallant/Desktop/for_github/ragtag_2022/pk_0.2_to_bbrach_final/intervals.list

docker exec gatk ./gatk SplitIntervals -R /jrgdata/ragtag.scaffold.fasta -L /jrgdata/intervals.list -scatter 50 -O /jrgdata/interval-files -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW

docker exec gatk cat /jrgdata/ragtag.scaffold.fasta.fai  | cut -f1,2 > /Users/jasongallant/Desktop/for_github/ragtag_2022/pk_0.2_to_bbrach_final/intervals.list

docker exec gatk ./gatk PreprocessIntervals -R /jrgdata/ragtag.scaffold.fasta --bin-length 200000 --padding 0 -O /jrgdata/preprocessed_intervals.interval_list

docker exec gatk ./gatk SplitIntervals -R /jrgdata/ragtag.scaffold.fasta -L /jrgdata/preprocessed_intervals.interval_list --scatter-count 100 -O /jrgdata/unpadded-interval-files --subdivision-mode INTERVAL_SUBDIVISION

cd ~/Desktop/for_github/pkings_firecloud_3/

cp -r ~/Desktop/for_github/ragtag_2022/pk_0.2_to_bbrach_final/interval-files .
cp -r ~/Desktop/for_github/ragtag_2022/pk_0.2_to_bbrach_final/unpadded-interval-files .

ls -1 ./interval-files > ./ifiles.txt

ls -1 ./unpadded-interval-files > ./uifiles.txt

awk '{ print "gs://fc-secure-001a4fd5-e796-46ba-8aa8-de834384e0bd/interval-files/" $0}' ./ifiles.txt > ./scattered-intervals-file.txt

awk '{ print "gs://fc-secure-001a4fd5-e796-46ba-8aa8-de834384e0bd/unpadded-interval-files/" $0}' ./uifiles.txt >./unpadded_intervals_file.txt

rm ./ifiles.txt
rm ./uifiles.txt

#### Modifications for GenomicsDB (comes later)

For some reason that I've never been able to figure out, the last 1-2 unpadded-interval-files cram many hundreds or thousands of contigs, depending on the contiguity of the assembly.  This causes the genotyping steps to fail later on because of an OOM issue.  The easiest way to fix this is to manually subdivided any unpadded intervals files that have >200 contigs in them to additional files.  If you forget to do this, you can do it after everything chokes without losing progress.  This can be done by doing the following:

cat 0098-scattered.intervals | grep -v "@" | split --lines 200
cat 0099-scattered.intervals.old | grep -v "@" | split --lines 200

You'll have to manually cut and paste the headers into these files and give them corresponding interval names > the last interval successfully ran in order to keep your progress.  Upload these, and modify the "unpadded_intervals_file.txt" appropriately, and everything should work fine


#### Upload to Workspace Bucket
1. Transfer `./unpadded_intervals_file.txt`, `./scattered-intervals-file.txt` and the directories `interval-files` and `unpadded-interval-files` to the Google Bucket

## Setup Workspace

1. Upload `./data_model/participant_2022.txt`to Terra
2. Upload `./data_model/sample_2022.tsv`, to Terra, which is generated from Airtable Inventory, modified in excel to reflect paths to reads uploaded in previous section
3. Upload `./data_model/pkings_gwas_3-workspace-attributes.tsv` to Terra
4. Run Workflow `./workflows/01_FastqToSam_v1.wdl` with configuration `./configurations/01_FastqToSam.1.json` to convert Fastq Files to uBAMs

## Merge uBAMs
1. Download The "Samples.TSV" from FireCloud
2. Run this code:
    awk -F "\t" '{print $13 >> ("./sample_metadata/unaligned_bams_" $12 ".txt")} ; close("./sample_metadata/unaligned_bams_" $12".txt")' ~/Downloads/Samples.TSV
3. Upload `sample_metadata` directory to bucket
4. Download Participants.TSV file from Terra and modify with column "unaligned_bam_list" pointing to google bucket location generated in step 3 above
5. Run `./workflows/collect_unaligned_bams_by_participant.1.wdl` with configuration `./configurations/collect_unaligned_bams_by_participant.json`

## Run Haplotype Caller
1. Run `./workflows/haplotype_caller.2.wdl` with configuration `./configurations/haplotype_caller.2.json`

This step takes a long time and is probably the most expensive part of the pipeline.  Estimate for 177 samples is ~24h

## Run Joint Genotyping
1. Run `./workflows/msu_efishlab_joint-discovery-gatk4.3.wdl` with configuration `msu_efishlab_joint-discovery-gatk4.3.json`
This step takes a long time Estimate for 177 samples is ~12-24h

You now have a RAW VCF for downstream analysis!!
