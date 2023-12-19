## Copyright Broad Institute, 2018
## Modified by JRG for Terra Nov 2019
##
## This WDL implements the joint discovery and VQSR filtering portion of the GATK
## Best Practices (June 2016) for germline SNP and Indel discovery in human
## whole-genome sequencing (WGS) and exome sequencing data.
##
## Requirements/expectations :
## - One or more GVCFs produced by HaplotypeCaller in GVCF mode
## - Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
##
## Outputs :
## - A VCF file and its index, filtered using variant quality score recalibration
##   (VQSR) with genotypes for all samples present in the input VCF. All sites that
##   are present in the input VCF are retained; filtered sites are annotated as such
##   in the FILTER field.
## - Note that the sample_names is what the sample will be called in the output, but not necessarily what the sample name is called in its GVCF.
##
##
## Cromwell version support
## - Successfully tested on v31
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow JointGenotyping {
  # Input Sample
  String callset_name

  Array[String] sample_names
  Array[File] input_gvcfs
  Array[File] input_gvcfs_indices

  # Reference and Resources
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File jg_mod_unpadded_intervals_file

  # Runtime attributes
  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.3.0.0"])
  String? gatk_path_override
  String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])

  Int? small_disk_override
  Int small_disk = select_first([small_disk_override, "100"])
  Int? medium_disk_override
  Int medium_disk = select_first([medium_disk_override, "200"])
  Int? large_disk_override
  Int large_disk = select_first([large_disk_override, "300"])
  Int? huge_disk_override
  Int huge_disk = select_first([huge_disk_override, "400"])

  String? preemptible_tries_override
  Int preemptible_tries = select_first([preemptible_tries_override, "3"])

  Array[String] unpadded_intervals = read_lines(jg_mod_unpadded_intervals_file)
  
  scatter (idx in range(length(unpadded_intervals))) {
    # the batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call ImportGVCFs {
      input:
        sample_names = sample_names,
        input_gvcfs = input_gvcfs,
        input_gvcfs_indices = input_gvcfs_indices,
        interval = unpadded_intervals[idx],
        workspace_dir_name = "genomicsdb",
        disk_size = huge_disk,
        batch_size = 50,
        docker = gatk_docker,
        gatk_path = gatk_path,
        preemptible_tries = preemptible_tries
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        disk_size = huge_disk,
        docker = gatk_docker,
        gatk_path = gatk_path,
        preemptible_tries = preemptible_tries
    }

  }

# for small callsets we can gather the VCF shards and then collect metrics on it
  call GatherVcfs as FinalGatherVcf {
    input:
      input_vcfs_fofn = write_lines(GenotypeGVCFs.output_vcf),
      output_vcf_name = callset_name + ".raw.vcf.gz",
      disk_size = huge_disk,
      docker = gatk_docker,
      gatk_path = gatk_path,
      preemptible_tries = preemptible_tries
  }


  output {
    # outputs from the small callset path through the wdl
    File? output_vcf = FinalGatherVcf.output_vcf
    File? output_vcf_index = FinalGatherVcf.output_vcf_index

    # output the interval list generated/used by this run workflow
  }
}

task GetNumberOfSamples {
  File sample_name_map
  String docker
  Int preemptible_tries
  command <<<
    wc -l ${sample_name_map} | awk '{print $1}'
  >>>
  runtime {
    docker: docker
    memory: "1 GB"
    preemptible: preemptible_tries
  }
  output {
    Int sample_count = read_int(stdout())
  }
}

task ImportGVCFs {
  Array[String] sample_names
  Array[String] input_gvcfs
  Array[String] input_gvcfs_indices
  String interval

  String workspace_dir_name

  String gatk_path
  String docker
  Int disk_size
  Int preemptible_tries
  Int batch_size

  command <<<
    set -e
    set -o pipefail

    python << CODE
    gvcfs = ['${sep="','" input_gvcfs}']
    sample_names = ['${sep="','" sample_names}']

    if len(gvcfs)!= len(sample_names):
      exit(1)

    with open("inputs.list", "w") as fi:
      for i in range(len(gvcfs)):
        fi.write(sample_names[i] + "\t" + gvcfs[i] + "\n")

    CODE

    rm -rf ${workspace_dir_name}

    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    
    if [[ ${interval} == *"remaining-interval-files"* ]]; then

        ${gatk_path} --java-options "-Xms8000m -Xmx25000m" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ${workspace_dir_name} \
        --batch-size ${batch_size} \
        -L ${interval} \
        --sample-name-map inputs.list \
        --reader-threads 5 \
        --merge-input-intervals \
        --consolidate \
        --merge-contigs-into-num-partitions 1

    else

        ${gatk_path} --java-options "-Xms8000m -Xmx25000m" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ${workspace_dir_name} \
        --batch-size ${batch_size} \
        -L ${interval} \
        --sample-name-map inputs.list \
        --reader-threads 5 \
        --merge-input-intervals \
        --bypass-feature-reader \
        --consolidate
        
    fi

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

  >>>
  runtime {
    docker: docker
    memory: "30 GB"
    cpu: "4"
    disks: "local-disk " + disk_size + " SSD"
    preemptible: preemptible_tries
  }
  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File workspace_tar
  String interval

  String output_vcf_filename

  String gatk_path

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String docker
  Int disk_size
  Int preemptible_tries

  command <<<
    set -e

    tar -xf ${workspace_tar}
    WORKSPACE=$( basename ${workspace_tar} .tar)

    ${gatk_path} --java-options "-Xms8000m -Xmx25000m" \
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --use-new-qual-calculator \
     -V gendb://$WORKSPACE \
     -all-sites \
     -L ${interval}
  >>>
  runtime {
    docker: docker
    memory: "60 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " SSD"
    preemptible: preemptible_tries
  }
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename
  String gatk_path

  String docker
  Int disk_size
  Int preemptible_tries

  command {
    set -e

    ${gatk_path} --java-options "-Xmx3g -Xms3g" \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

  }
  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
  }
}

task GatherVcfs {
  File input_vcfs_fofn
  String output_vcf_name
  String gatk_path

  String docker
  Int disk_size
  Int preemptible_tries

  command <<<
    set -e

    # Now using NIO to localize the vcfs but the input file must have a ".list" extension
    mv ${input_vcfs_fofn} inputs.list

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    ${gatk_path} --java-options "-Xmx6g -Xms6g" \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    --input inputs.list \
    --output ${output_vcf_name}

    ${gatk_path} --java-options "-Xmx6g -Xms6g" \
    IndexFeatureFile \
    -I ${output_vcf_name}
  >>>
  runtime {
    docker: docker
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + 1000 + " SSD"
    preemptible: preemptible_tries
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task DynamicallyCombineIntervals {
  File intervals
  Int merge_count
  Int preemptible_tries

  command {
    python << CODE
    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

    count = 0
    chain_count = ${merge_count}
    l_chr, l_start, l_end = "", 0, 0
    lines_to_write = []
    with open("${intervals}") as f:
        with open("out.intervals", "w") as f1:
            for line in f.readlines():
                # initialization
                if count == 0:
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue
                # reached number to combine, so spit out and start over
                if count == chain_count:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue

                c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                if c_chr == w_chr and c_start == w_end + 1:
                    w_end = c_end
                    count += 1
                    continue
                # not adjacent, end here and start a new chain
                else:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
            if l_char != w_chr or l_start != w_start or l_end != w_end:
                add_interval(w_chr, w_start, w_end)
            f1.writelines("\n".join(lines_to_write))
    CODE
  }

  runtime {
    memory: "3 GB"
    preemptible: preemptible_tries
    docker: "python:2.7"
  }

  output {
    File output_intervals = "out.intervals"
  }
}
