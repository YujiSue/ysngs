version 1.0
# FastQC
task fastqc {
    input {
        String fq
        String dir
    }
    command <<<
        mkdir -p ~{dir}
        $HYM_APP/FastQC/fastqc -o ~{dir} ~{fq} &
    >>>
    output {}
}

task cutadapt1 {
    input {
        File fq
        String dir
    }
    command <<<
        cutadapt -b ATCACCGACTGCCCATAGAGAGGCTGAGAC ~{fq} > $HYM_DATA/~{dir}/$HYM_DATA/~{dir}/sub(basename(~{fq}), ".fq", ".cut.fq")
        echo $HYM_DATA/~{dir}/sub(basename(~{fq}), ".fq", ".cut.fq")
    >>>
    output {
        String cut = read_string(stdout())
    }
}
task cutadapt2 {
    input {
        Array[String] fq
        String dir
    }
    command <<<

    >>>
    output {
        String cut = read_string(stdout())
    }
}
# fastp
task fastp {
    input {
        Boolean paired = false
        Array[String] fq
        String in1 = "-i ~{fq[0]}"
        String in2 = if paired then "-I ~{fq[1]}" else ""
        String dir
        String out
        String out1 = if paired then "-o ~{dir}/~{out}_1.fq" else "-o ~{dir}/~{out}.fq"
        String out2 = if paired then "-O ~{dir}/~{out}_2.fq" else ""
        String? option = ""
        Int? thread = 2
    }
    command <<<
        mkdir -p $HYM_DATA/~{dir}
        $HYM_APP/fastp -w ~{thread} \
          ~{option} \
          ~{in1} \
          ~{in2} \
          ~{out1} \
          ~{out2} \
          -h ~{dir}/~{out}.html \
          -j ~{dir}/~{out}.json
        ls ~{dir}/~{out}*.fq
    >>>
    output {
        Array[String] filtered = read_lines(stdout())
    }
}
task fastp1 {
    input {
        Array[String] fq
        String dir
        String ext
        Map[String, String] condition        
    }
    command <<<
        fastp -i ~{fq[0]} \
            -o ~{dir}/sub(basename(~{fq[0]}), ~{ext}, ".filtered.fq") \
            -l ~{condition["min-len"]} \
            -e ~{condition["min-qual"]} \
            -f ~{condition["trim-head"]} \
            -b ~{condition["trim-tail"]} \
            -a ~{condition["cut-adapt"]} \
            -h $HYM_DATA/~{dir}/sub(basename(~{fq[0]}), ~{ext}, ".html")
        echo $HYM_DATA/~{dir}/sub(basename(~{fq[0]}), ~{ext}, ".filtered.fq")
    >>>
    output {
        Array[String] qcfq = glob(read_string(stdout()))
    }
}
task fastp2 {
    input {
        Array[String] fq
        String dir
        String ext
        Map[String, String] condition        
    }
    command <<<
        
    >>>
    output {}
}
task filteredfq {
    input {
        String dir
    }
    command <<<
        echo $HYM_DATA/~{dir}/*.filtered.fq
    >>>
    output {
        Array[String] filtered = glob(read_string(stdout()))
    }
}