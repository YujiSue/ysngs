version 1.0
import "common.wdl" as common
import "samtools.wdl" as samtools
import "bwa.wdl" as bwa
import "bowtie2.wdl" as bowtie2
import "gatk.wdl" as gatk
import "star.wdl" as star
import "hisat.wdl" as hisat
import "rsem.wdl" as rsem
workflow fa2idx {
    input {
        Boolean useBWA
        Boolean useBowtie2
        Boolean useGATK
        Boolean useSTAR
        Boolean useHISAT2
        Boolean useRSEM
        String fa
        String label
        Int thread
    }
    call samtools.makefaidx {
        input:
            fa = fa
    }
    if (useBWA) {
        call bwa.bwaindex {
            input:
                fa = fa,
                label = label
        }
    }
    if (useBowtie2) {
        call bowtie2.bowtindex {
            input:
                fa = fa,
                label = label,
                thread = thread
        }
    }
    if (useGATK) {
        call gatk.refindex {
            input:
                fa = fa,
                label = label
        }
    }
    
}