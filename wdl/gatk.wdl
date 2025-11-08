version 1.0

# Make reference index
task gatkindex {
  input {
    String fa
    String refdir = sub(fa, "/[^/]+$", "")
    String label
  }
  command <<<
    $HYM_APP/gatk/gatk CreateSequenceDictionary -R ~{fa} -O ~{refdir}/~{label}.dict
  >>>
  output {}
}
# Make intervals
task gatkinterval {
  input {
    String fa
    String refdir = sub(fa, "/[^/]+$", "")
    String interval
  }
  command <<<
    $HYM_APP/gatk/gatk PreprocessIntervals \
      -R ~{fa} \
      --padding 0 \
      -imr OVERLAPPING_ONLY \
      -O ~{refdir}/~{interval}
  >>>
  output {}
}
# Mark duplicate
task markdp {
  input {
    String bam
    String dir
    String name
    String out = "~{dir}/~{name}.bam"
  }
  command <<< 
    java -jar $HYM_APP/picard.jar MarkDuplicates \
      -I ~{bam} \
      -O ~{out} \
      -M ~{dir}/~{name}.metric.txt
  >>>
  output {
    String aligned = out
  }
}
# Base recalibration
task brecal {
  input {
    String reference
    String bam
    Array[String] sites

    String dir
    String name
    String out = "~{dir}/~{name}.bq.table"

    Int ram = 1
  }
  command <<<
    mkdir -p ~{dir}
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" BaseRecalibrator \
      -R ~{reference} \
      -I ~{bam} \
      --known-sites ~{sep=" --known-sites " sites} \
      -O ~{out}
  >>>
  output {
    String bq = out
  }
}
# Apply base recalibration
task applybq {
  input {
    String reference
    String bam
    String table
    
    String dir
    String name
    String out = "~{dir}/~{name}.recal.bam"

    Int ram = 1    
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx{ram}g" ApplyBQSR \
      -R ~{reference} \
      -I ~{bam} \
      --bqsr-recal-file ~{table} \
      -O ~{out}
  >>>
  output {
    String recalibrated = out
  }
}
# Haplotype call
task haplocall {
  input {
    String reference
    String bam

    Boolean target = false
    String bed = ""
    String region = if target then "-L ~{bed}" else ""

    String dir
    String name
    String out = "~{dir}/~{name}.gatk.g.vcf"

    Int thread = 2
    Int ram = 1
  }
  command <<<
    mkdir -p ~{dir}
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" HaplotypeCaller \
      -R ~{reference} \
      ~{region} \
      -I ~{bam} \
      -O ~{out} \
      -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
      --native-pair-hmm-threads ~{thread}
  >>>
  output {
    String total = out
  }
}
# Genotype
task genotype {
  input {
    String reference
    String gvcf

    String dir
    String name

    String out = "~{dir}/~{name}.vcf"

    Int ram = 1
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" GenotypeGVCFs \
      -R ~{reference} \
      -V ~{gvcf} \
      -O ~{out}
  >>>
  output {
    String variants = out
  }
}
# Variant quality recalbration
task vrecal {
  input {
    String reference

    Array[String] filters
    String resources
    
    String mode
    
    String vcf

    String dir
    String name
    
    Int ram = 1
  }
  command <<<
    mkdir -p ~{dir}/others
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" VariantRecalibrator \
      -R ~{reference} \
      -V ~{vcf} \
      ~{resources} \
      ~{sep=" -an " filters} \
      -mode ~{mode} \
      -O ~{dir}/others/~{name}.recal \
      --tranches-file ~{dir}/others/~{name}.tranches \
      --rscript-file ~{dir}/others/~{name}.plot.R
    echo ~{dir}/others/~{name}
  >>>
  output {
    String vqpre = read_string(stdout())
  }
}
# Apply base recalibration
task applyvq {
  input {
    String reference

    String prefix
    String mode
    String level

    String vcf

    String dir
    String name
    String out = "~{dir}/~{name}.~{mode}.recal.vcf"
    
    Int ram = 1
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx{ram}g" ApplyVQSR \
      -R ~{reference} \
      -V ~{vcf} \
      --recal-file ~{prefix}.recal \
      --tranches-file ~{prefix}.tranches \
      --truth-sensitivity-filter-level ~{level} \
      -mode ~{mode} \
      -O ~{out}
  >>>
  output {
    String recalibrated = out
  }
}
# Variant filter
task gatkfilter {
  input {
    String reference
    String mode

    String vcf
    String condition

    String dir
    String name
    String out = "~{dir}/~{name}.~{mode}.filtered.vcf"

    Int ram = 1
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" VariantFiltration \
      -R ~{reference} \
      -V ~{vcf} \
      -O ~{out} \
      ~{condition}
  >>>
  output {
    String filtered = out
  }
}
# Split variant list
task vcfsplit {
  input {
    String reference
    String vcf

    String dir
    String name
    String snp_vcf_path = "~{dir}/variant/~{name}.snp.vcf"
    String indel_vcf_path = "~{dir}/variant/~{name}.indel.vcf"

    Int ram = 1
  }
  
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" SelectVariants \
      -R ~{reference} \
      -V ~{vcf} \
      --select-type-to-include SNP \
      -O ~{snp_vcf_path}
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" SelectVariants \
      -R ~{reference} \
      -V ~{vcf} \
      --select-type-to-include INDEL \
      -O ~{indel_vcf_path}
  >>>
  output {
    String snp_vcf = snp_vcf_path
    String indel_vcf = indel_vcf_path
  }
}
# Concat variant list
task vcfconcat {
  input {
    Array[String] vcfs
    
    String dir
    String name
    String out = "~{dir}/~{name}.merged.vcf"

    Int ram = 1
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" MergeVcfs \
      -I ~{sep=" -I " vcfs} \
      -O ~{out}
  >>>
  output {
    String merged = out
  }
}

# VC workflow
workflow gatkvc {
    input {
        String reference
        String bam
        
        Boolean target_seq = false
        String target_region = ""

        String out_dir
        String out_name
        
        Boolean base_recal = false
        Array[String] known_sites = []

        Boolean var_recal = false
        Array[String] snp_recal_filter = []
        Array[Pair[String,String]] snp_sources = [] 
        Map[String,String] snp_sources_opts = {}
        Array[String] indel_recal_filter = []
        Array[Pair[String,String]] indel_sources = [] 
        Map[String,String] indel_sources_opts = {}

        Boolean var_filter = false
        Array[Pair[String,String]] snp_filter = []
        Array[Pair[String,String]] indel_filter = [] 
        
        Int thread = 2
        Int ram = 1
    }
    # Base score recalibration
    if (base_recal) {
        # Make recal. table
        call brecal {
            input:
                reference = reference,
                bam = bam,
                dir = out_dir + "/others",
                name = out_name,
                sites = known_sites,
                ram = ram
        }
        # Apply recal. table
        call applybq {
            input:
                reference = reference,
                bam = bam,
                table = brecal.bq,
                dir = out_dir + "/align",
                name = out_name,
                ram = ram
        }
    }
    # Reset input BAM
    String bam2 = select_first([applybq.recalibrated, bam])
    
    # Haplotype caller
    call haplocall {
        input:
            reference = reference,
            bam = bam2,
            target = target_seq,
            bed = target_region,
            dir = out_dir + "/variant",
            name = out_name,
            thread = thread,
            ram = ram
    }
    # 
    call genotype {
        input:
            ram = ram,
            reference = reference,
            gvcf = haplocall.total,
            dir = out_dir + "/variant",
            name = out_name
    }
    # Variant qual. recalibration
    if (var_recal) {
        ## For SNP
        ### Interpret resource list
        scatter (r in snp_sources) {
          String snp_src = "--resource:" + r.left + "," + snp_sources_opts[r.left] + " " + r.right
        }
        String snp_src_opt = "~{sep=' ' snp_src}"
        ### Training
        call vrecal as snprecal {
            input:
                ram = ram,
                reference = reference,
                filters = snp_recal_filter,
                resources = snp_src_opt,
                mode = 'SNP',
                vcf = genotype.variants,
                dir = out_dir,
                name = out_name
        }
        ### Apply
        call applyvq as snpvq {
            input: 
                ram = ram,
                reference = reference,
                prefix = snprecal.vqpre,
                mode = 'SNP',
                vcf = genotype.variants,
                dir = out_dir,
                name = out_name,
                level = '99.0'
        }
        ## For INDEL
        ### Interpret resource list
        scatter (r in indel_sources) {
          String indel_src = "--resource:" + r.left + "," + indel_sources_opts[r.left] + " " + r.right
        }
        String indel_src_opt = "~{sep=' ' indel_src}"
        ### Training
        call vrecal as indelrecal {
            input:
                ram = ram,
                reference = reference,
                filters = indel_recal_filter,
                resources = indel_src_opt,
                mode = 'INDEL',
                vcf = genotype.variants,
                dir = out_dir,
                name = out_name
        }        
        ### Apply
        call applyvq as indelvq {
            input: 
                ram = ram,
                reference = reference,
                prefix = indelrecal.vqpre,
                mode = 'INDEL',
                vcf = snpvq.recalibrated,
                dir = out_dir,
                name = out_name,
                level = '99.0'
        }
    }
    # Hard filter
    if (var_filter) {
        # Split into SNP/INDEL list 
        call vcfsplit {
            input:
                ram = ram,
                reference = reference,
                vcf = genotype.variants,
                dir = out_dir,
                name = out_name
        }
        ## Make SNP filter condition
        scatter (f in snp_filter) {
          String snp_filt = "--filter-expression \"" + f.right + "\" --filter-name \"" + f.left + "\""
        }
        String snp_filt_opt = "~{sep=' ' snp_filt}"
        ## Filter
        call gatkfilter as snpfilter {
            input:
                ram = ram,
                reference = reference,
                vcf = vcfsplit.snp_vcf,
                mode = "snp",
                dir = out_dir,
                name = out_name,
                condition = snp_filt_opt
        }
        ## Make INDEL filter condition
        scatter (f in indel_filter) {
          String indel_filt = "--filter-expression \"" + f.right + "\" --filter-name \"" + f.left + "\""
        }
        String indel_filt_opt = "~{sep=' ' indel_filt}"
        
        ## Filter
        call gatkfilter as indelfilter {
            input:
                ram = ram,
                reference = reference,
                vcf = vcfsplit.indel_vcf,
                mode = "indel",
                dir = out_dir,
                name = out_name,
                condition = indel_filt_opt
        }
        # Merge filtered variants
        call vcfconcat {
            input:
                ram = ram,
                vcfs = [snpfilter.filtered, indelfilter.filtered],
                dir = out_dir,
                name = out_name
        }
    }
    #
    String out = select_first([vcfconcat.merged, indelvq.recalibrated, genotype.variants])
    output {
        String gvcf = haplocall.total
        String vcf = out
    }
}
