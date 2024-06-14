version 1.0
import "samtools.wdl"
workflow bam2fq {
    input {
        String bam
        String out_dir
    }
    call samtools.bam2fq {
        input: 
          bam = bam,
          dir = out_dir
    }
    output {
        String fq = bam2fq.fq
    } 
}
