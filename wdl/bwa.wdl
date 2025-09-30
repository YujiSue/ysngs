version 1.0
task bwamap {
    input {
        Array[String] fq
        String info
        String ref
        String dir
        String name
        Int thread
    }
    command <<< 
        $HYM_APP/bwa mem -t ~{thread} -Y -M \
          -R "@RG\t~{info}" \
          ~{ref} ~{sep=" " fq} > ~{dir}/~{name}.bwa.sam
        echo ~{dir}/~{name}.bwa.sam
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
    >>>
    output {}
}
