version 1.0
# Mapping
task bowtmap {
  input {
    Array[String] fq
    Boolean paired
    String f1 = if paired then "-1 ~{fq[0]}" else "-U ~{fq[0]}"
    String f2 = if paired then "-2 ~{fq[1]}" else ""

    String smplid
    String sample
    String library
    String platform
    
    String ref
    
    String dir
    String name
    String out = "~{dir}/~{name}.bowt.sam"

    Int thread
  }
  command <<< 
    mkdir -p ~{dir}
    $HYM_APP/bowtie2/bowtie2 \
      -p ~{thread} \
      --rg-id ~{smplid} \
      --rg SM:~{sample} \
      --rg LB:~{library} \
      --rg PL:~{platform} \
      ~{f1} \
      ~{f2} \
      -x ~{ref} \
      -S ~{out}
  >>>
  output {
    String sam = out
  }
}
# Make index
task bowtindex {
    input {
        String fa
        String label
        Int thread = 2
    }
    command <<< 
        cd $HYM_REF
        $HYM_APP/bowtie2/bowtie2-build --threads ~{thread} -f ~{fa} ~{label}
        cd $HYM_WS
    >>>
    output {}
}