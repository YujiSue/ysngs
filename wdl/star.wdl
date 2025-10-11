version 1.0
# mapping
task starmap {
  input {
    Array[String] fq
    String ref
    
    String dir
    String name
    String out = "~{dir}/~{name}.bam"
    
    Int thread
  }
  command <<< 
    $HYM_APP/STAR/bin/Linux_x86_64/STAR \
      --genomeDir ~{ref} \
      --outFileNamePrefix ~{dir}/~{name} \
      --quantMode TranscriptomeSAM \
      --outSAMtype BAM SortedByCoordinate \
      --readFilesIn ~{sep=" " fq} \
      --runThreadN ~{thread}
  >>>
  output {
    String aligned = out
  }
}
# Make index
task starindex {
  input {
    String fa
    String gtf
	  String dir
    Int thread
  }
  command <<< 
    mkdir -p ~{dir}
    $HYM_APP/STAR \
      --runMode genomeGenerate \
      --genomeDir ~{dir} \
      --genomeFastaFiles ~{fa} \
      --sjdbGTFfile ~{gtf} \
      --runThreadN ~{thread}
  >>>
  output {}
}
