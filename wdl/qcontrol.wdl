version 1.0
import "common.wdl" as common
import "qc.wdl" as qc

workflow qcontrol {
    input {
        Boolean paired
        Array[String] fastq
        String out_name
        String out_dir
    }
    # Primary check
    scatter (fq in fastq) {
      call qc.fastqc as primary_check {
        input:
          fq = fq,
          dir = "~{out_dir}/QC/primary"
      }
    }
    # Cut adapter and filter reads
    call qc.fastp as fastp { 
        input: 
            paired = paired,
            fq = fastq,
            dir = "~{out_dir}/QC",
            out = out_name
    }
    # Final check
    scatter (fq in fastp.filtered) {
      call qc.fastqc as final_check {
        input:
          fq = fq,
          dir = "~{out_dir}/QC/final"
      }
    }
    output {
        Array[String] filtered = fastp.filtered
    }
}