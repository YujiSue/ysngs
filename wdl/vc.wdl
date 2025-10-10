version 1.0
import "common.wdl" as common
import "samtools.wdl" as samtools
import "gatk.wdl" as gatk
import "google.wdl" as google

workflow varcall {
    input {
        String species = ""
        String reference
        Boolean target_seq = false
        String target_region = ""
        
        Array[String] known_sites = []
        Array[String] snp_recal_filter = ["QD", "MQ", "MQRankSum", "ReadPosRankSum", "FS", "SOR"]
        Array[Pair[String,String]] snp_sources = [] 
        Map[String,String] snp_sources_opts = {}
        Array[String] indel_recal_filter = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR"]
        Array[Pair[String,String]] indel_sources = [] 
        Map[String,String] indel_sources_opts = {}

        String bam

        String out_dir
        String out_name

        Boolean use_gatk = true
        Boolean base_recal = false
        Boolean var_recal = false
        Boolean var_filter = false

        Boolean use_bcftools = false

        Boolean use_gdv = false
        Boolean gdv_gpu = false
        String gdv_ver = ""
        String gdv_image = ""

        Array[Pair[String,String]] snp_filter = [
            ("LowQD", "QD < 2.0"),
            ("HighFS", "FS > 60.0"),
            ("HighSOR", "SOR > 3.0"),
            ("LowMQ", "MQ < 40.0"),
            ("LowMQRankSum", "MQRankSum < -12.5"),
            ("LowReadPosRankSum", "ReadPosRankSum < -8.0")
        ]
        Array[Pair[String,String]] indel_filter = [
            ("LowQD", "QD < 2.0"),
            ("HighFS", "FS > 200.0"),
            ("LowReadPosRankSum", "ReadPosRankSum < -20.0")
        ]

        Int thread
        Int ram
    }
    # GATK VC
    if (use_gatk) {
        call gatk.gatkvc {
            input:
                reference = reference,
                bam = bam,
                target_seq = target_seq,
                target_region = target_region,
                
                out_dir = out_dir,
                out_name = out_name,
                
                base_recal = base_recal,
                known_sites = known_sites,

                var_recal = var_recal,
                snp_recal_filter = snp_recal_filter,
                snp_sources = snp_sources,
                snp_sources_opts = snp_sources_opts,
                
                indel_recal_filter = indel_recal_filter,
                indel_sources = indel_sources,
                indel_sources_opts = indel_sources_opts,
                
                var_filter = var_filter,
                snp_filter = snp_filter,
                indel_filter = indel_filter,
                
                thread = thread,
                ram = ram
        }
    }
    # BCFtools mpileup
    if (use_bcftools) {
        call samtools.bcfcall {
            input:
                reference = reference,
                bam = bam,
                dir = out_dir,
                name = out_name
        }
    }
    # Google Deepvariant
    if (use_gdv) {
        call google.gdvcall {
            input:
                use_gpu = gdv_gpu,
                ver = gdv_ver,
                path = gdv_image,
                reference = reference,
                target_seq = target_seq,
                target = target_region,
                bam = bam,
                dir = out_dir,
                name = out_name,
                thread = thread
        }        
    }
    String pre = select_first([gatkvc.gvcf, gdvcall.gvcf, ""])
    String out = select_first([gatkvc.vcf, bcfcall.vcf, gdvcall.vcf])
    output {
        String gvcf = pre
        String vcf = out
    }
}
