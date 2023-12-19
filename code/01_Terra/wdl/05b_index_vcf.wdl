workflow IndexVCF {

  String output_vcf

  # Runtime attributes
  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.3.0.0"])
  String? gatk_path_override
  String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])

  Int? huge_disk_override
  Int huge_disk = select_first([huge_disk_override, "400"])

  String? preemptible_tries_override
  Int preemptible_tries = select_first([preemptible_tries_override, "3"])
  
  
# for small callsets we can gather the VCF shards and then collect metrics on it
  call IndexVCFs as IndexFinalVCF {
    input:
      output_vcf_name = output_vcf,
      disk_size = huge_disk,
      docker = gatk_docker,
      gatk_path = gatk_path,
      preemptible_tries = preemptible_tries
  }


  output {
    # outputs from the small callset path through the wdl
    File? output_vcf_index = IndexFinalVCF.output_vcf_index

    # output the interval list generated/used by this run workflow
  }
}

task IndexVCFs {
  String output_vcf_name
  String gatk_path

  String docker
  Int disk_size
  Int preemptible_tries

  command <<<
    set -e

    ${gatk_path} --java-options "-Xmx6g -Xms6g" \
    IndexFeatureFile \
    --I ${output_vcf_name}
  >>>
  runtime {
    docker: docker
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " SSD"
    preemptible: preemptible_tries
  }
  output {
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}