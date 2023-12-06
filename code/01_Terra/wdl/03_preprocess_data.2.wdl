## Copyright Broad Institute, 2019
## 
## This WDL pipeline implements data pre-processing according to the GATK Best Practices.  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
## Software version requirements 
## - GATK 4 or later
## - BWA 0.7.15-r1140
## - Picard 2.16.0-SNAPSHOT
## - Samtools 1.3.1 (using htslib 1.3.1)
## - Python 2.7
##
## Cromwell version support 
## - Successfully tested on v37
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow PreProcessingForVariantDiscovery_GATK4 {

 String sample_name
 String ref_name

 File flowcell_unmapped_bams_list
 String unmapped_bam_suffix
 
 File ref_fasta
 File ref_fasta_index
 File ref_dict
 
 String? bwa_commandline_override
 String bwa_commandline = select_first([bwa_commandline_override, "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"]) 
 Int compression_level
 
 String? gatk_docker_override
 String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.1.0.0"])
 String? gatk_path_override
 String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])

 String? gotc_docker_override
 String gotc_docker = select_first([gotc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"])
 String? gotc_path_override
 String gotc_path = select_first([gotc_path_override, "/usr/gitc/"])

 String? python_docker_override
 String python_docker = select_first([python_docker_override, "python:2.7"])  

 Int flowcell_small_disk
 Int flowcell_medium_disk
 Int agg_small_disk
 Int agg_medium_disk
 Int agg_large_disk

 String? preemptible_tries_override
 Int preemptible_tries = select_first([preemptible_tries_override, "3"])

 String base_file_name = sample_name + "." + ref_name

 Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)

 # Get the version of BWA to include in the PG record in the header of the BAM produced 
 # by MergeBamAlignment. 
 call GetBwaVersion {
   input: 
     docker_image = gotc_docker,
     bwa_path = gotc_path,
     preemptible_tries = preemptible_tries
 }

 # Align flowcell-level unmapped input bams in parallel
 scatter (unmapped_bam in flowcell_unmapped_bams) {

   # Get the basename, i.e. strip the filepath and the extension
   String bam_basename = basename(unmapped_bam, unmapped_bam_suffix)

   # Map reads to reference
   call SamToFastqAndBwaMem {
     input:
       input_bam = unmapped_bam,
       bwa_commandline = bwa_commandline,
       output_bam_basename = bam_basename + ".unmerged",
       ref_fasta = ref_fasta,
       ref_fasta_index = ref_fasta_index,
       ref_dict = ref_dict,
       docker_image = gotc_docker,
       bwa_path = gotc_path,
       gotc_path = gotc_path,
       disk_size = flowcell_medium_disk,
       preemptible_tries = preemptible_tries,
       compression_level = compression_level
    }

   # Merge original uBAM and BWA-aligned BAM 
   call MergeBamAlignment {
     input:
       unmapped_bam = unmapped_bam,
       bwa_commandline = bwa_commandline,
       bwa_version = GetBwaVersion.version,
       aligned_bam = SamToFastqAndBwaMem.output_bam,
       output_bam_basename = bam_basename + ".aligned.unsorted",
       ref_fasta = ref_fasta,
       ref_fasta_index = ref_fasta_index,
       ref_dict = ref_dict,
       docker_image = gatk_docker,
       gatk_path = gatk_path,
       disk_size = flowcell_medium_disk,
       preemptible_tries = preemptible_tries,
       compression_level = compression_level
   }
 }

 # Aggregate aligned+merged flowcell BAM files and mark duplicates
 # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
 # to avoid having to spend time just merging BAM files.
 call MarkDuplicates {
   input:
     input_bams = MergeBamAlignment.output_bam,
     output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
     metrics_filename = base_file_name + ".duplicate_metrics",
     docker_image = gatk_docker,
     gatk_path = gatk_path,
     disk_size = agg_large_disk,
     compression_level = compression_level,
     preemptible_tries = preemptible_tries
 }

 # Sort aggregated+deduped BAM file and fix tags
 call SortAndFixTags {
   input:
     input_bam = MarkDuplicates.output_bam,
     output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
     ref_dict = ref_dict,
     ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index,
     docker_image = gatk_docker,
     gatk_path = gatk_path,
     disk_size = agg_large_disk,
     preemptible_tries = 0,
     compression_level = compression_level
 }

 # Outputs that will be retained when execution is complete  
 output {
   File duplication_metrics = MarkDuplicates.duplicate_metrics
   File analysis_ready_bam = SortAndFixTags.output_bam
   File analysis_ready_bam_index = SortAndFixTags.output_bam_index
   File analysis_ready_bam_md5 = SortAndFixTags.output_bam_md5
 } 
}

# TASK DEFINITIONS

# Get version of BWA
task GetBwaVersion {
 
 Int preemptible_tries
 String mem_size

 String docker_image
 String bwa_path

 command {
   # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
   # because the sed may also fail with that error and that is something we actually want to fail on.
   ${bwa_path}bwa 2>&1 | \
   grep -e '^Version' | \
   sed 's/Version: //'
 }
 runtime {
   preemptible: preemptible_tries
   docker: docker_image
   memory: mem_size
 }
 output {
   String version = read_string(stdout())
 }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
 File input_bam
 String bwa_commandline
 String output_bam_basename
 File ref_fasta
 File ref_fasta_index
 File ref_dict

 # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
 # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
 # references such as b37 and hg19.
 File? ref_alt
 File ref_amb
 File ref_ann
 File ref_bwt
 File ref_pac
 File ref_sa

 Int compression_level
 Int preemptible_tries
 Int disk_size
 String mem_size
 String num_cpu

 String docker_image
 String bwa_path
 String gotc_path
 String java_opt

 command <<<
   set -o pipefail
   set -e

   # set the bash variable needed for the command-line
   bash_ref_fasta=${ref_fasta}

   java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${gotc_path}picard.jar \
     SamToFastq \
     INPUT=${input_bam} \
     FASTQ=/dev/stdout \
     INTERLEAVE=true \
     NON_PF=true \
   | \
   ${bwa_path}${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
   | \
   samtools view -1 - > ${output_bam_basename}.bam

 >>>
 runtime {
   preemptible: preemptible_tries
   docker: docker_image
   memory: mem_size
   cpu: num_cpu
   disks: "local-disk " + disk_size + " HDD"
 }
 output {
   File output_bam = "${output_bam_basename}.bam"
   File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
 }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
 File unmapped_bam
 String bwa_commandline
 String bwa_version
 File aligned_bam
 String output_bam_basename
 File ref_fasta
 File ref_fasta_index
 File ref_dict

 Int compression_level
 Int preemptible_tries
 Int disk_size
 String mem_size

 String docker_image
 String gatk_path
 String java_opt

 command {
   # set the bash variable needed for the command-line
   bash_ref_fasta=${ref_fasta}
   ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
     MergeBamAlignment \
     --VALIDATION_STRINGENCY SILENT \
     --EXPECTED_ORIENTATIONS FR \
     --ATTRIBUTES_TO_RETAIN X0 \
     --ALIGNED_BAM ${aligned_bam} \
     --UNMAPPED_BAM ${unmapped_bam} \
     --OUTPUT ${output_bam_basename}.bam \
     --REFERENCE_SEQUENCE ${ref_fasta} \
     --PAIRED_RUN true \
     --SORT_ORDER "unsorted" \
     --IS_BISULFITE_SEQUENCE false \
     --ALIGNED_READS_ONLY false \
     --CLIP_ADAPTERS false \
     --MAX_RECORDS_IN_RAM 2000000 \
     --ADD_MATE_CIGAR true \
     --MAX_INSERTIONS_OR_DELETIONS -1 \
     --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
     --PROGRAM_RECORD_ID "bwamem" \
     --PROGRAM_GROUP_VERSION "${bwa_version}" \
     --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
     --PROGRAM_GROUP_NAME "bwamem" \
     --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
     --ALIGNER_PROPER_PAIR_FLAGS true \
     --UNMAP_CONTAMINANT_READS true
 }
 runtime {
   preemptible: preemptible_tries
   docker: docker_image
   memory: mem_size
   disks: "local-disk " + disk_size + " HDD"
 }
 output {
   File output_bam = "${output_bam_basename}.bam"
 }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
 File input_bam
 String output_bam_basename
 File ref_dict
 File ref_fasta
 File ref_fasta_index
 
 Int compression_level
 Int preemptible_tries
 Int disk_size
 String mem_size

 String docker_image
 String gatk_path
 String java_opt_sort
 String java_opt_fix

 command {
   set -o pipefail

   ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_sort}" \
     SortSam \
     --INPUT ${input_bam} \
     --OUTPUT /dev/stdout \
     --SORT_ORDER "coordinate" \
     --CREATE_INDEX false \
     --CREATE_MD5_FILE false \
   | \
   ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_fix}" \
     SetNmMdAndUqTags \
     --INPUT /dev/stdin \
     --OUTPUT ${output_bam_basename}.bam \
     --CREATE_INDEX true \
     --CREATE_MD5_FILE true \
     --REFERENCE_SEQUENCE ${ref_fasta}
 }
 runtime {
   preemptible: preemptible_tries
   docker: docker_image
   memory: mem_size
   disks: "local-disk " + disk_size + " HDD"
 }
 output {
   File output_bam = "${output_bam_basename}.bam"
   File output_bam_index = "${output_bam_basename}.bai"
   File output_bam_md5 = "${output_bam_basename}.bam.md5"
 }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
 Array[File] input_bams
 String output_bam_basename
 String metrics_filename
 
 Int compression_level
 Int preemptible_tries
 Int disk_size
 String mem_size

 String docker_image
 String gatk_path
 String java_opt

# Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
# This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
# While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
 command {
   ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
     MarkDuplicates \
     --INPUT ${sep=' --INPUT ' input_bams} \
     --OUTPUT ${output_bam_basename}.bam \
     --METRICS_FILE ${metrics_filename} \
     --VALIDATION_STRINGENCY SILENT \
     --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
     --ASSUME_SORT_ORDER "queryname" \
     --CREATE_MD5_FILE true
 }
 runtime {
   preemptible: preemptible_tries
   docker: docker_image
   memory: mem_size
   disks: "local-disk " + disk_size + " HDD"
 }
 output {
   File output_bam = "${output_bam_basename}.bam"
   File duplicate_metrics = "${metrics_filename}"
 }
}