version 1.0
import "common.wdl" as common
import "samtools.wdl" as samtools
import "bwa.wdl" as bwa
import "gatk.wdl" as gatk
workflow bwaalign {
    input {
        Array[String] fastq
        String read_info
        String reference
        String out_dir
        String out_name
        Int thread = 2
        Boolean clean = true
    }
    call bwa.bwamap {
        input:
            fq = fastq,
            info = read_info,
            ref = reference,
            dir =  out_dir,
            name = out_name,
            thread = thread
    }
    call samtools.sam2bam {
        input:
            sam = bwamap.sam,
            name = out_name,
            dir = out_dir,
            thread = thread
    }
    call samtools.bamsort {
        input:
            bam = sam2bam.rawbam,
            name = out_name,
            dir = out_dir,
            thread = thread
    }
    call gatk.markdp {
        input:
            bam = bamsort.sorted,
            name = out_name,
            dir = out_dir,
            name = out_name
    }
    call samtools.makeindex {
        input:
            bam = markdp.aligned
    }
    if (clean) {
        call common.remove {
            input:
                src = [
                    bwamap.sam,
                    sam2bam.rawbam,
                    bamsort.sorted
                ]
        }
    }
    output {
        String aligned = markdp.aligned
    }
}
