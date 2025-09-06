version 1.0
task bwamap {
    input {
        Array[String] fq
        Map[String, String] read_info
        String dir
        String ref
        Int thread
    }
    command <<< 
        $HYM_APP/bwa mem -t ~{thread} -Y -M \
          -R "@RG\tID:~{read_info['id']}\tSM:~{read_info['sample']}\tPL:~{read_info['platform']}" \
          ~{ref} ~{sep=" " fq} > ~{dir}/aligned.bwa.sam
        echo ~{dir}/aligned.bwa.sam
    >>>
    output {
        String sam = read_string(stdout())
    }
}
task bwaindex {
    input {
        String fa
        String label
    }
    command <<< 
        $HYM_APP/bwa index -p ~{label} ~{fa}
#        mv ~{label}.amb $HYM_REF
#        mv ~{label}.ann $HYM_REF
#        mv ~{label}.bwt $HYM_REF
#        mv ~{label}.pac $HYM_REF
#        mv ~{label}.sa $HYM_REF
    >>>
    output {}
}
