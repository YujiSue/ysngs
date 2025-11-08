version 1.0

import "common.wdl" as common
import "qc.wdl" as qc

workflow qcontrol {
    input {
        Array[String] fastq
        Boolean paired

        String out_name
        String out_dir

        Int thread = 2
    }
    # Primary check
    scatter (fq in fastq) {
      call qc.fastqc as primary_qc {
        input:
          fq = fq,
          dir = out_dir + "/QC/primary"
      }
    }
    # Cut adapter and filter reads
    call qc.fastp { 
        input: 
            paired = paired,
            fq = fastq,
            dir = out_dir + "/QC",
            name = out_name,
            thread = thread
    }
    # Final check
    scatter (fq in fastp.filtered) {
      call qc.fastqc as secondary_qc {
        input:
          fq = fq,
          dir = out_dir + "/QC/final"
      }
    }
    output {
        Array[String] filtered = fastp.filtered
    }
}