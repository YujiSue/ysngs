version 1.0
task refindex {
  input {
    String fa
    String label
  }
  command <<<
    $HYM_APP/gatk/gatk CreateSequenceDictionary -R $HYM_REF/~{fa} -O $HYM_REF/~{label} 
  >>>
  output {}
}
task refinterval {
  input {
    File fa
    String interval
  }
  command <<<
    $HYM_APP/gatk/gatk PreprocessIntervals -R ~{fa} --padding 0 -imr OVERLAPPING_ONLY -O ~{interval}
  >>>
  output {}
}
task markdp {
  input {
    File bam
    String dir
    String name
  }
  command <<< 
    java -jar $HYM_APP/picard.jar MarkDuplicates -I ~{bam} -O $HYM_DATA/~{dir}/~{name}.bam -M $HYM_DATA/~{dir}/~{name}.metric.txt
    echo $HYM_DATA/~{dir}/~{name}.bam
  >>>
  output {
    String aligned = read_string(stdout())
  }
}
task haplocall {
  input {
    String bam_dir
    String bam
    String fa
    String out_dir
    Int ram
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" HaplotypeCaller \
      -R $HYM_REF/~{fa} \
      -I $HYM_DATA/~{bam_dir}/~{bam} \
      -O $HYM_DATA/~{out_dir}/sub(basename(bam), ".bam", ".gatk.g.vcf") \
      -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
    echo $HYM_DATA/~{out_dir}/sub(basename(bam), ".bam", ".gatk.g.vcf")
  >>>
  output {
    String gvcf = read_string(stdout())
  }
}
task haplocallt {
  input {
    String bam_dir
    String bam
    String fa
    String out_dir
    File bed
    Int ram
  }
  command <<<
    $HYM_APP/gatk/gatk --java-options "-Xmx~{ram}g" HaplotypeCaller \
      -R $HYM_REF/~{fa} \
      -L ~{bed} \
      -I $HYM_DATA/~{bam_dir}/~{bam} \
      -O $HYM_DATA/~{out_dir}/sub(basename(bam), ".bam", ".gatk.g.vcf") \
      -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
    echo $HYM_DATA/~{out_dir}/sub(basename(bam), ".bam", ".gatk.g.vcf")
  >>>
  output {
    String gvcf = read_string(stdout())
  }  
}
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

