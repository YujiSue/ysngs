version 1.0
import "common.wdl" as common
import "samtools.wdl" as samtools
import "bwa.wdl" as bwa
import "gatk.wdl" as gatk
workflow bwaalign {
    input {
        Array[String] fastq
        Map[String, String] read_info
        String ref_label
        String out_dir
        String out_name
        Int thread
        Boolean clean
    }
    call bwa.bwamap {
        input:
            fq = fastq,
            read_info = read_info,
            ref = ref_label,
            thread = thread
    }
    call samtools.sam2bam {
        input:
            sam = bwamap.sam,
            thread = thread
    }
    call samtools.bamsort {
        input:
            bam = sam2bam.rawbam,
            thread = thread
    }
    call gatk.markdp {
        input:
            bam = bamsort.sorted,
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
