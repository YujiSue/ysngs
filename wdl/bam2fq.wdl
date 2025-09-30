version 1.0
import "samtools.wdl"
workflow bam2fq {
    input {
        String bam
        String out_dir
        Boolean paired
    }
    if (paired) {
        call samtools.bam2fq {
            input: 
                bam = bam,
                dir = out_dir
        }
    } 
    if (!paired) {
        call samtools.bam2pair {
            input: 
                bam = bam,
                dir = out_dir
        }
    }
    output {} 
}
