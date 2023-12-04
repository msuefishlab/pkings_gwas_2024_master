workflow generateunalignedBAM {
  call FastqtoSam
}

task FastqtoSam {
  File R1
  File R2
  String sampleName
  String LB
  String Lane
  String Barcode
  String SC
  String PL
    
  Int memoryGb
  Int diskSpaceGb

  command {
java -Dsamjdk.buffer_size=131072 \
-Dsamjdk.use_async_io=true \
-Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8G \
-jar /usr/gitc/picard.jar FastqToSam \
FASTQ=${R1} \
FASTQ2=${R2} \
OUTPUT=${sampleName}_${LB}.${Lane}.${Barcode}_fastqtosam.bam \
READ_GROUP_NAME=${LB}.${Lane} \
SAMPLE_NAME=${sampleName} \
LIBRARY_NAME=${LB} \
PLATFORM_UNIT=${LB}.${Lane} \
PLATFORM=${PL} \
SEQUENCING_CENTER=${SC}
    }

  output {
  File unaligned_bam = "${sampleName}_${LB}.${Lane}.${Barcode}_fastqtosam.bam"
  }

  runtime {
  docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
  memory: "${memoryGb} GB"
  cpu: "1"
  disks: "local-disk ${diskSpaceGb} HDD"
  
  }
}