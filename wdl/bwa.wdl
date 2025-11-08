version 1.0

task bwamap {
    input {
        Array[String] fq
        
        Boolean rg = true
        String? smplid = ""
        String? sample = ""
        String? library = ""
        String? platform = ""
        String readinfo = if rg then "-R \"@RG\\tID:~{smplid}\\tSM:~{sample}\\tLB:~{library}\\tPL:~{platform}\"" else ""

        String ref
        String dir
        String name
        String out = "~{dir}/~{name}.bwa.sam"
        Int? thread = 2
    }
    command <<< 
        mkdir -p ~{dir}
        $HYM_APP/bwa mem -t ~{thread} -Y -M \
          ~{readinfo} \
          ~{ref} \
          ~{sep=" " fq} > ~{out}
    >>>
    output {
        String sam = out
    }
}
task bwaindex {
    input {
        String fa
        String dir
        String label
    }
    command <<< 
        $HYM_APP/bwa index -p ~{dir}/~{label} ~{fa}
    >>>
    output {}
}
