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
# Cutadapt (Old code. Use fastp)
task cutadapt {
    input {
        Array[String] fq
        Boolean paired

        Boolean cut3end = true
        String adapter3
        String cut3 = if cut3end then (if paired then "-a ~{adapter3} -A ~{adapter3}" else "-a ~{adapter3}") else ""

        Boolean cut5end = false
        String adapter5 = ""
        String cut5 = if cut5end then (if paired then "-g ~{adapter5} -G ~{adapter5}" else "-g ~{adapter5}") else ""

        
        String dir
        String name

        String out1 = if paired then "-o ~{dir}/~{name}.cut_1.fq" else "-o ~{dir}/~{name}.cut.fq"
        String out2 = if paired then "-p ~{dir}/~{name}.cut_2.fq" else ""
    }
    command <<<
        cutadapt \
          ~{cut3} \
          ~{cut5} \
          ~{out1} \
          ~{out2} \
          ~{sep=" " fq}          
        ls ~{dir}/~{name}.cut*.fq
    >>>
    output {
        Array[String] cut = read_lines(stdout())
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
        String name
        String f1 = if paired then "~{dir}/~{name}_1.fq" else "~{dir}/~{name}.fq"
        String f2 = if paired then "~{dir}/~{name}_2.fq" else ""
        String out1 = "-o ~{f1}"
        String out2 = "-O ~{f2}"
        
        
        String option = ""
        Int thread = 2
    }
    command <<<
        mkdir -p ~{dir}
        $HYM_APP/fastp -w ~{thread} \
          ~{option} \
          ~{in1} \
          ~{in2} \
          ~{out1} \
          ~{out2} \
          -h ~{dir}/~{name}.html \
          -j ~{dir}/~{name}.json
    >>>
    output {
        Array[String] filtered = if paired then [f1] else [f1, f2]
    }
}