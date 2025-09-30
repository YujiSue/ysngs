version 1.0
import "common.wdl" as common
import "gatk.wdl" as gatk

workflow gatk_vc {
    input {
        Boolean? base_recal = false
        Boolean? var_recal = false
        Boolean? var_filter = false
        
        Array[String] bams
        String reference
        String out_name
        String out_dir
    }
    scatter (bam in bams) {
        String in = bam
        if (base_recal) {
            call gatk.brecal as brecal {}
            call gatk.applybq as applybq {}
            
            in = applybq.recalibrated
        }
        call gatk.haplocall as haplocall {
            input:
                bam = in,

        }
        call gatk.genotype as genotype {

        }
        if (var_recal) {
            call gatk.brecal as brecal {}
            call gatk.applybq as applybq {}
            
            in = applybq.recalibrated
        }
        if (var_filter) {
            #call gatk.
        }

    }
    output {
        Array[String] variants = fastp.filtered
    }
}

workflow gatk_vc {

    
}