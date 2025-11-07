version 1.0

task hisatmap {
    input {
        Array[String] fq
        Boolean paired
        String f1 = if paired then "-1 ~{fq[0]}" else "-U ~{fq[0]}"
        String f2 = if paired then "-2 ~{fq[1]}" else ""

        String ref
        
        String dir
        String name
        String out = "~{dir}/~{name}.hisat.sam"

        Int thread
    }
    command <<<
        $HYM_APP/hisat2/hisat2 \
          -x ~{ref} \
          ~{f1} \
          ~{f2} \          
          -p ~{thread} \
          -S ~{out}
    >>>
    output {
        String sam = out
    }
}
# Make suppl. data
task hisatexon {
    input {
        String gtf
        String dir
        String name
        String out = "~{dir}/~{name}.exons.txt"
    }
    command <<< 
        mkdir -p ~{dir}/~{name}
        python $HYM_APP/hisat2/scripts/hisat2_extract_exons.py ~{gtf} > ~{out}
        
    >>>
    output {
        String exons = out
    }
}
task hisatss {
    input {
        String gtf
        String dir
        String name
        String out = "~{dir}/~{name}.ss.txt"
    }
    command <<< 
        mkdir -p ~{dir}/~{name}
        python $HYM_APP/hisat2/scripts/hisat2_extract_splice_sites.py ~{gtf} > ~{out}        
    >>>
    output {
        String ss = out
    }
}
# make index
task hisatindex {
    input {
		Boolean rnaseq = false
        String fa
        String gtf = ''
        String exons = ''
        String ss = ''
        String option = if rnaseq then "--exons ~{exons} --ss ~{ss}" else ""

        String dir
        String label
        Int thread
    }
    command <<< 
        mkdir -p ~{dir}/~{label}
        $HYM_APP/hisat2/hisat2-build \
          -p ~{thread} \
          ~{option}
          ~{fa} \
          ~{dir}/~{label}
    >>>
    output {}
}