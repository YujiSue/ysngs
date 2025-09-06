version 1.0
import "common.wdl" as common
import "qc.wdl" as qc

workflow qcheck {
    input {
        Array[String] fastq
        String out_dir
    }
    scatter (file in fastq) {
        call qc.fastqc { 
            input: 
                fq = file,
                dir = out_dir
        }
    }
}