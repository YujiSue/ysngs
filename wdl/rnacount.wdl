version 1.0
# feature count
task fcount {
    input {
        String bam
        Boolean paired = false
        String opt = if paired then "-p -B -C" else ""

        String gtf
        String type = "exon"
        String topt = if type == "" then "" else "-t ~{type}"
        String by = "gene_id"
        String bopt = if by == "" then "" else "-g ~{by}"

        String dir
        String name
        String out = "~{dir}/count/~{name}.fcount.txt"
        
        Int thread = 2
    }
    command <<<
        featureCounts \
          ~{opt} \
          ~{topt} \
          ~{bopt} \
          -T ~{thread} \
          -a ~{gtf} \
          -o ~{out} \
          ~{bam}
    >>>
    output {
        String count = out
    }
}

# HTSeq
task htscount {
    input {
        String bam
        Boolean paired = false
        
        String gtf
        String type = "exon"
        String topt = if type == "" then "" else "-t ~{type}"
        String by = "gene_id"
        String bopt = if by == "" then "" else "-i ~{by}"

        String dir
        String name
        String out = "~{dir}/count/~{name}.htseq.txt"
    }
    command <<<
        htseq-count \
          -f bam \
          -r name \
          -s no \
          ~{topt} \
          ~{bopt} \
          ~{bam} \
          ~{gtf} > ~{out}
    >>>
    output {
        String count = out
    }
}

# stringtie
## Find novel isoform
task stiefind {
    input {
        String bam
        String gtf
        String dir
        String name
        String out = "~{dir}/~{name}.gtf"
        Int thread
    }
    command <<<
        stringtie \
          ~{bam} \
          -p ~{thread} \
          -G ~{gtf} \
          -o ~{out}
    >>>
    output {
        String novel = out
    }
}
## Make new GTF to count novel isoform
task stiemerge {
    input {
        String gtf
        String list
        String dir
        String out = "~{dir}/stie.merged.gtf"
    }
    command <<<
        stringtie --merge \
          -G ~{gtf} \
          -o ~{out} \
          ~{list}
    >>>
    output {
        String merged = out
    }
}
## Count
task stiecount {
    input {
        String bam
        String gtf
        String dir
        String name
        String out = "~{dir}/stringtie/~{name}"
        
        Int thread = 2
    }
    command <<<
        mkdir ~{out}
        stringtie ~{bam} \
          -G ~{gtf} \
          -o ~{out}/output.gtf \
          -p ~{thread}
    >>>
    output {
        String result = out
    }
}
## Conversion
task stieexport {
    input {
        String list
        String dir
        String name
        String out1 = "~{dir}/count/{name}.stie.gene.count.csv"
        String out2 = "~{dir}/count/{name}.stie.transcripts.count.csv"
    }
    command <<<
        python $HYM_APP/stringtie/prepDE.py \
          -i ~{list} \
          -g ~{out1} \
          -t ~{out2}
    >>>
    output {
        String count = out1
        String count2 = out2
    }
}
