version 1.0
import "common.wdl" as common
import "cellranger.wdl" as cr

workflow sctranscriptome {
    input {
        Boolean multi

        Boolean export_bam = true
        String bam = if export_bam then "true" else "false"

        String reference = ""
        String config = ""
        String dir = ""
	String subdir = ""
        String file_prefix = ""
        String out_name = ""
        
        Boolean use_scanpy = false
        Boolean use_seurat = false
        
        #Boolean novel_transcript = false

        Int thread = 0
        Int ram = 0
    }
    # Counting
    if (!multi) {
        call cr.sccount {
            input:
                reference = reference,
                bam = bam,
                dir = dir,
		subdir = subdir,
                prefix = file_prefix,
                outname = out_name,
                core = thread,
                ram = ram
        }
    }
    if (multi) {
        call cr.scmulti {
            input:
                config = config,
		dir = dir,
                outname = out_name,
                core = thread,
                ram = ram
        }
    }
    String matrix = select_first([sccount.result, scmulti.result, ""])
    if (use_scanpy) {
        
    }
    if (use_seurat) {
        
    }
    output {}
}

