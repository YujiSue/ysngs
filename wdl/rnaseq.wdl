version 1.0
import "common.wdl" as common

import "qcontrol.wdl" as qc
import "mapping.wdl" as map
#import "cuff.wdl" as cuff
import "rnacount.wdl" as counter
import "rsem.wdl" as rsem

workflow transcriptome {
    input {
        Boolean single_cell = false
        
        Boolean mapping = true
        Boolean use_bowtie = false
        Boolean use_star = false
        Boolean use_hisat = false
        Boolean rsem_map = false
        Boolean export_bam = true

        Array[String] fastq = []
        Boolean paired = false
        String map_reference = ""

        String? bam
        String out_dir
        String out_name
        
        Boolean use_fcount = false
        Boolean use_htseq = false
        Boolean use_rsem = false
        String ref_fa = ""
        String ref_gtf = ""
        
        Boolean novel_isoform = false
        #Boolean use_cuff = false
        Boolean use_stie = false

        Int thread = 2
    }
    # QC
    call qc.qcontrol {
        input:
            paired = paired,
            fastq = fastq,
            out_dir = out_dir,
            out_name = out_name,
            thread = thread
    }
    # Mapping
    if (mapping) {
        if (rsem_map) {
            call rsem.rsemcount {
                input:
                    mapping = true,
                    export_bam = export_bam,
                    use_bowtie = use_bowtie,
                    use_star = use_star,
                    use_hisat = use_hisat,
                    paired = paired,
                    fq = qcontrol.filtered,
                    dir = out_dir,
                    name = out_name,
                    ref = map_reference,
                    thread = thread
            }
        }
        if (!rsem_map) {
            if (use_bowtie) {}
            if (use_star) {}
            if (use_hisat) {}            
        }
        
    }
    String genome_bam = select_first([rsemcount.gbam,bam])
    String transcript_bam = select_first([rsemcount.tbam,bam])
    if (use_fcount) {
        call counter.fcount {
            input:
                bam = genome_bam,
                paired = paired,
                gtf = ref_gtf,
                dir = out_dir,
                name = out_name,
                thread = thread
        }
    }
    if (use_htseq) {
        call counter.htscount {
            input:
                bam = genome_bam,
                paired = paired,
                gtf = ref_gtf,
                dir = out_dir,
                name = out_name
        }
    }
    output {}
}

