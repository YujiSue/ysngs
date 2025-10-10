version 1.0

import "common.wdl" as common
import "samtools.wdl" as samtools
import "bwa.wdl" as bwa
import "bowtie2.wdl" as bowt
import "gatk.wdl" as gatk

workflow align {
    input {
        Boolean use_bwa = false
        Boolean use_bowtie = false
        Boolean clean = true
        
        Array[String] fastq
        Map[String, String] read_info
        Boolean paired

        String reference

        String out_dir
        String out_name

        Int thread = 2
    }
    # BWA
    if (use_bwa) {
        call bwa.bwamap {
            input:
                fq = fastq,
                
                smplid = read_info["id"],
                sample = read_info["sample"],
                library = read_info["lib"],
                platform = read_info["platform"],

                ref = reference,
                dir =  out_dir + "/align",
                name = out_name,

                thread = thread
        }
    }
    # Bowtie2
    if (use_bowtie) {
        call bowt.bowtmap {
            input:
                fq = fastq,
                paired = paired,

                smplid = read_info["id"],
                sample = read_info["sample"],
                library = read_info["lib"],
                platform = read_info["platform"],

                ref = reference,
                dir =  out_dir + "/align",
                name = out_name,

                thread = thread
        }
    }
    #
    call samtools.sam2bam {
        input:
            sam = select_first([bwamap.sam, bowtmap.sam]),
            name = out_name+".raw",
            dir = out_dir + "/align",
            thread = thread
    }
    call samtools.bamsort {
        input:
            bam = sam2bam.rawbam,
            name = out_name+".sorted",
            dir = out_dir + "/align",
            thread = thread
    }
    call gatk.markdp {
        input:
            bam = bamsort.sorted,
            name = out_name,
            dir = out_dir + "/align",
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
                    select_first([bwamap.sam, bowtmap.sam]),
                    sam2bam.rawbam,
                    bamsort.sorted
                ],
                exclude = [
                    markdp.aligned,
                    makeindex.index
                ]
        }
    }
    output {
        String aligned = markdp.aligned
    }
}