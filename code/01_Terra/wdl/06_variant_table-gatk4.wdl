workflow MakeVariantTable {

  String input_vcf
  Array[String] var_types

  # Reference and Resources
  String ref_fasta

  # Runtime attributes
  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.3.0.0"])
  String? gatk_path_override
  String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])

  Int? huge_disk_override
  Int huge_disk = select_first([huge_disk_override, "400"])

  String? preemptible_tries_override
  Int preemptible_tries = select_first([preemptible_tries_override, "3"])
  
  
  call GetVariantTable  {
    input:
      input_vcf = input_vcf,
      ref_fasta = ref_fasta,
      disk_size = huge_disk,
      docker = gatk_docker,
      gatk_path = gatk_path,
      preemptible_tries = preemptible_tries
  }
  scatter (vtype in var_types) {
    call SplitVarType{
      input:
        vtype = vtype,
        table = GetVariantTable.vartable,
        docker = gatk_docker,
        disk_size = huge_disk,
        preemptible_tries = preemptible_tries
    }
  }

  output {
    Array[File] splitTables = SplitVarType.splitvtable
  }
}

task GetVariantTable {
  String input_vcf
  String gatk_path
  
  String ref_fasta

  String docker
  Int disk_size
  Int preemptible_tries

  command <<<
    set -e

    gatk VariantsToTable \
    -R ${ref_fasta} \
    -V ${input_vcf} \
    -F CHROM \
    -F POS \
    -F TYPE \
    -F ID \
    -F QUAL \
    -F AC \
    -F HET \
    -F HOM-REF \
    -F HOM-VAR \
    -F AN \
    -F BaseQRankSum \
    -F DP \
    -F FS \
    -F MQ \
    -F MQRankSum \
    -F QD \
    -F ReadPosRankSum \
    -F SOR \
    -GF GQ \
    -O ${input_vcf}.tab
    
  >>>
  runtime {
    docker: docker
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " SSD"
    preemptible: preemptible_tries
  }
  output {
    File vartable = "${input_vcf}.tab"
  }
}

task SplitVarType {
  
  String vtype
  String table

  String docker
  Int disk_size
  Int preemptible_tries
  
  command {
    awk 'NR<=3 || $3 == "'${vtype}'"' ${table} > ${table}.${vtype}.tab
  }
  
  runtime {
    docker: docker
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " SSD"
    preemptible: preemptible_tries
  }

  output{
    File splitvtable="${table}.${vtype}.tab"
  }
}
