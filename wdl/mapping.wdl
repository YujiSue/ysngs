version 1.0

import "common.wdl" as common
import "samtools.wdl" as samtools
import "bwa.wdl" as bwa
import "bowtie2.wdl" as bowt
import "star.wdl" as star
import "hisat.wdl" as hisat
import "gatk.wdl" as gatk

workflow align {
    input {
        String mode = "genome"

        Boolean use_bwa = false
        Boolean use_bowtie = false
        Boolean use_star = false
        Boolean use_hisat = false
        Boolean clean = true
        
        Array[String] fastq
        Map[String, String] read_info
        Boolean paired

        String reference

        String out_dir
        String out_name

        Int thread = 2
    }
    if (mode == "genome") {
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
        # SAM=>BAM
        call samtools.sam2bam {
            input:
                sam = select_first([bwamap.sam, bowtmap.sam]),
                name = out_name+".raw",
                dir = out_dir + "/align",
                thread = thread
        }
        # Sort
        call samtools.bamsort {
            input:
                bam = sam2bam.rawbam,
                name = out_name+".sorted",
                dir = out_dir + "/align",
                thread = thread
        }
        # Detect duplication
        call gatk.markdp {
            input:
                bam = bamsort.sorted,
                name = out_name,
                dir = out_dir + "/align",
                name = out_name
        }
        # Index
        call samtools.makeindex {
            input:
                bam = markdp.aligned
        }
        # Clear intermediates
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
    }
    if (mode == "rna") {
        if (use_bowtie) {
            
        }
        if (use_star) {
            call star.starmap {
                input:
                    fq = fastq,
                    ref = reference,
                    dir = out_dir + "/align",
                    name = out_name,
                    thread = thread
            }
        }
        if (use_hisat) {
            
        }
    }
    output {
        String aligned = if mode == "genome" then markdp.aligned else select_first([starmap.aligned])
    }
}