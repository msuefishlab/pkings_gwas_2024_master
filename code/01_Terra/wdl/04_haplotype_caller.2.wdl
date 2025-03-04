## Copyright Broad Institute, 2017
##
## This WDL workflow runs HaplotypeCaller from GATK4 in GVCF mode on a single sample
## according to the GATK Best Practices (June 2016), scattered across intervals.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##
## Cromwell version support
## - Successfully tested on v29
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
workflow HaplotypeCallerGvcf_GATK4 {
  File input_bam
  File input_bam_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File intervals_list

  String picard_docker
  String gatk_docker

  String picard_path
  String gatk_path

  Array[String] intervals_array = read_lines(intervals_list)

  String sample_basename = basename(input_bam, ".bam")

  String gvcf_name = sample_basename + ".g.vcf.gz"
  String gvcf_index = sample_basename + ".g.vcf.gz.tbi"

  # Call variants in parallel over grouped calling intervals
  scatter (the_interval in intervals_array) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        interval = the_interval,
        gvcf_name = gvcf_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }
  }
  #Array[File] gvcf_intervals = HaplotypeCaller.output_gvcf

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      vcf_name = gvcf_name,
      vcf_index = gvcf_index,
      docker_image = picard_docker,
      picard_path = picard_path
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_merged_gvcf = MergeGVCFs.output_vcf
    File output_merged_gvcf_index = MergeGVCFs.output_vcf_index
  }
}



# TASK DEFINITIONS

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  String gvcf_name
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String interval
  Int? interval_padding
  Float? contamination
  Int? max_alt_alleles

  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {
    ${gatk_path} --java-options ${java_opt} \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${gvcf_name} \
      -L ${interval} \
      -ip ${default=100 interval_padding} \
      -contamination ${default=0 contamination} \
      --max-alternate-alleles ${default=3 max_alt_alleles} \
      -ERC GVCF
  }

  runtime {
    docker: docker_image
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_gvcf = "${gvcf_name}"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  Array [File] input_vcfs
  String vcf_name
  String vcf_index

  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String picard_path
  String java_opt

  command {
    java ${java_opt} -jar ${picard_path}picard.jar \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${vcf_name}
  }

  runtime {
    docker: docker_image
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
}

  output {
    File output_vcf = "${vcf_name}"
    File output_vcf_index = "${vcf_index}"
  }
}
