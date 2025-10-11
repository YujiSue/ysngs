version 1.0
import "common.wdl" as common
import "samtools.wdl" as samtools
import "bwa.wdl" as bwa
import "bowtie2.wdl" as bowtie2
import "gatk.wdl" as gatk
import "star.wdl" as star
import "hisat.wdl" as hisat
import "rsem.wdl" as rsem
# 
workflow fa2idx {
    input {
        Boolean use_bwa = false
        Boolean use_bowtie = false
        Boolean use_gatk = false
        Boolean use_star = false
        Boolean use_hisat = false
        Boolean use_rsem = false
        Boolean rsem_map = false
        
        String ref_fasta
        String ref_label
        String ref_gtf = ""
        
        Int thread  =2
    }
    call samtools.makefaidx {
        input:
            fa = ref_fasta
    }
    if (use_bwa) {
        call bwa.bwaindex {
            input:
                fa = ref_fasta,
                label = ref_label
        }
    }
    if (use_bowtie) {
        call bowtie2.bowtindex {
            input:
                fa = ref_fasta,
                label = ref_label,
                thread = thread
        }
    }
    if (use_gatk) {
        call gatk.gatkindex {
            input:
                fa = ref_fasta,
                label = ref_label
        }
    }
    if (use_star) {
        call star.starindex {
            input:
                fa = ref_fasta,
                gtf = ref_gtf,
                dir = out_dir,
                thread = thread
        }
    }
    if (use_hisat) {
        call hisat.hisatexon {
            input:
                gtf = ref_gtf,
                dir = out_dir,
                name = ref_label
        }
        call hisat.hisatss {
            input:
                gtf = ref_gtf,
                dir = out_dir,
                name = ref_label
        }
        call hisat.hisatindex {
            input:
                fa = ref_fasta,
                gtf = ref_gtf,
                exons = hisatexon.exons,
                ss = hisatss.ss,
                dir = out_dir,
                label = ref_label,
                thread = thread
        }
    }
}