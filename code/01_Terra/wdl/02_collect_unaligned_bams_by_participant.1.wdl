## Collect Unaligned BAMS by Participant
## Jason Gallant
##
## This WDL collects unaligned bam files by participant for downstream analysis.
## run on this.participants, should be set at the participant level

workflow unaliged_bam_array {
   File file_of_files
   Array[File] array_of_files = read_lines(file_of_files)

   output {
    Array[File] array_output = array_of_files
   }
}