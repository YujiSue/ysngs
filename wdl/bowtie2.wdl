version 1.0
task bowtmap {
  input {
    Array[File] fq
    Map[String, String] read_info
    String ref
    Int thread
  }
  command <<< 
    $HYM_APP/bowtie2/bowtie2 \
    --rg-id ~{read_info['id']} \
    --rg SM:~{read_info['sample']} \
    --rg PL:~{read_info['platform']} \
    -U ~{fq[0]} \
    -p ~{thread} \
    -x ~{ref} \
    -S aligned.bowt.sam
  >>>
  output {
    File sam = "aligned.bowt.sam"
  }
}
task bowtmap2 {
  input {
    Array[File] fq
    Map[String, String] read_info
    String ref
    Int thread
  }
  command <<< 
    bowtie2 --rg-id ~{read_info['id']} \
    --rg SM:~{read_info['sample']} \
    --rg PL:~{read_info['platform']} \
    -1 ~{fq[0]} -2 ~{fq[1]} \
    -p ~{thread} -x ~{ref} -S aligned.bowt.sam
  >>>
  output {
    File sam = "aligned.bowt.sam"
  }
}
task bowtindex {
    input {
        File fa
        String label
        Int thread
    }
    command <<< 
        $HYM_APP/bowtie2/bowtie2-build --threads ~{thread} -f ~{fa} $HYM_REF/~{label}
    >>>
    output {}
}
