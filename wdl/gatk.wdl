version 1.0
# Make reference index
task refindex {
  input {
    String fa
    String label
  }
  command <<<
    $HYM_APP/gatk/gatk CreateSequenceDictionary -R ~{fa} -O ~{label} 
  >>>
  output {}
}
# Make intervals
task refinterval {
  input {
    String fa
    String interval
  }
  command <<<
    $HYM_APP/gatk/gatk PreprocessIntervals -R ~{fa} --padding 0 -imr OVERLAPPING_ONLY -O ~{interval}
  >>>
  output {}
}
# Mark duplicate
task markdp {
  input {
    String bam
    String dir
    String name
  }
  command <<< 
    java -jar $HYM_APP/picard.jar MarkDuplicates -I ~{bam} -O ~{dir}/~{name}.bam -M ~{dir}/~{name}.metric.txt
    echo ~{dir}/~{name}.bam
  >>>
  output {
    String aligned = read_string(stdout())
  }
}
# Base recalibration
task brecal {
  input {
    Int? ram = 1
    String reference
    String dir
    String bam
    String name
    Array[String] sites
  }
  command <<<
    mkdir -p ~{dir}/others
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" BaseRecalibrator \
      -R ~{reference} \
      -I ~{bam} \
      ~{sep=" --known-sites " sites} \
      -O $HYM_DATA/~{dir}/others/~{name}.bq.table
    echo $HYM_DATA/~{dir}/others/~{name}.bq.table
  >>>
  output {
    String bq = read_string(stdout())
  }
}
# Apply base recalibration
task applybq {
  input {
    Int? ram = 1
    String reference
    String bam
    String table
    String dir
    String name
    Array[String] sites
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx{ram}g" ApplyBQSR \
      -R ~{reference} \
      -I ~{bam} \
      --bqsr-recal-file ~{table} \
      -O $HYM_DATA/~{dir}/aligned/~{name}.recal.bam
    echo $HYM_DATA/~{dir}/aligned/~{name}.recal.bam
  >>>
  output {
    String recalibrated = read_string(stdout())
  }
}
# Haplotype call
task haplocall {
  input {
    String bam_dir
    String bam
    String fa
    String? bed = ""
    String target = if bed == "" then "-L ~{bed}" else ""
    String out_dir
    Int? ram = 1
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" HaplotypeCaller \
      -R $HYM_REF/~{fa} \
      ~{target} \
      -I $HYM_DATA/~{bam_dir}/~{bam} \
      -O $HYM_DATA/~{out_dir}/sub(basename(bam), ".bam", ".gatk.g.vcf") \
      -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
    echo $HYM_DATA/~{out_dir}/sub(basename(bam), ".bam", ".gatk.g.vcf")
  >>>
  output {
    String gvcf = read_string(stdout())
  }
}
# Genotype
task genotype {
  input {
    String gvcf
    String out_dir
    String fa
    Int ram
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" GenotypeGVCFs \
      -R $HYM_REF/~{fa} -V ~{gvcf} -O $HYM_DATA/~{out_dir}/sub(basename(gvcf), ".g.vcf", ".vcf")
    echo $HYM_DATA/~{out_dir}/sub(basename(gvcf), ".g.vcf", ".vcf")
  >>>
  output {
    String vcf = read_string(stdout())
  }
}
# 
