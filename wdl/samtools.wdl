version 1.0
task bam2fq {
    input {
        String bam
        String name = sub(basename(bam), ".bam", ".fq") 
        String dir
    }
    command <<<
        mkdir -p ~{dir}
        samtools fastq ~{bam} > ~{dir}/~{name}
        echo ~{dir}/~{name}
    >>>
    output {
        String fq = read_string(stdout())
    }
}
task makefaidx {
    input {
        String fa
    }
    command <<<
        samtools faidx ~{fa}
        echo ~{fa}.fai
    >>>
    output {
        String fai = read_string(stdout())
    }
}
task sam2bam {
    input {
        String sam
        String dir
        Int thread
    }
    command <<<
        samtools view -@ ~{thread} \
          -b -o ~{dir}/raw.bam \
          ~{sam} 
        echo ~{dir}/raw.bam
    >>>
    output {
        String rawbam = read_string(stdout())
    }
}
task bamsort {
    input {
        String bam
        String dir
        String temp
        Int thread
    }
    command <<<
        samtools sort -l1 -T ~{temp} \
          -@ ~{thread} \
          -O bam \
          -o ~{dir}/sorted.bam \
          ~{bam}
        echo ~{dir}/sorted.bam
    >>>
    output {
        String sorted = read_string(stdout())
    }
}
task makeindex {
    input {
        String bam
    }
    command <<<
        samtools index ~{bam}
    >>>
    output {}
}