version 1.0

task starmap {
  input {
    Array[String] fq
    String out
    String ref
    Int thread
  }
  command <<< 
    $HYM_APP/STAR/bin/Linux_x86_64/STAR \
      --genomeDir $HYM_REF/STAR/~{ref} \
      --outFileNamePrefix $HYM_DATA/~{out}.star.bam \
      --quantMode TranscriptomeSAM \
      --outSAMtype BAM SortedByCoordinate
      --readFilesIn ~{sep=" " fq} \
      --runThreadN ~{thread}
    echo $HYM_DATA/~{out}.star.bam
  >>>
  output {
    String bam = read_string(stdout())
  }
}

task starindex {
  input {
    String fa
    String gtf
    String label
	Int thread
  }
  command <<< 
    mkdir -p $HYM_REF/STAR
    $HYM_APP/STAR/bin/Linux_x86_64/STAR \
      --runMode genomeGenerate \
      --genomeDir $HYM_REF/STAR/~{label} \
      --genomeFastaFiles $HYM_REF/~{fa} \
      --sjdbGTFfile $HYM_REF/~{gtf} \
      --runThreadN ~{thread}
    echo $HYM_REF/STAR/~{label}
  >>>
  output {
    String refdir = read_string(stdout())
  }
}
