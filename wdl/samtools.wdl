version 1.0
# BAM => single-read FASTQ
task bam2fq {
    input {
        String bam
        String name = basename(bam)
        String dir
    }
    command <<<
        mkdir -p ~{dir}
        samtools fastq ~{bam} > ~{dir}/~{name}.fq
        echo ~{dir}/~{name}
    >>>
    output {
        String fq = read_string(stdout())
    }
}
# BAM => paired-end FASTQ
task bam2pair {
    input {
        String bam
        String name = basename(bam)
        String dir
    }
    command <<<
        mkdir -p ~{dir}
        samtools collate -u -O ~{bam} | samtools fastq -1 ~{dir}/~{name}_1.fq -2 ~{dir}/~{name}_2.fq -0 /dev/null -s /dev/null -n
        ls ~{dir}/~{name}*.fq
    >>>
    output {
        Array[String] fq = read_lines(stdout())
    }
}
# Make fasta index
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
# SAM => BAM
task sam2bam {
    input {
        String sam
        String name = basename(sam)
        String dir
        Int? thread = 2
    }
    command <<<
        samtools view -@ ~{thread} \
          -b -o ~{dir}/~{name}.raw.bam \
          ~{sam} 
        echo ~{dir}/~{name}.raw.bam
    >>>
    output {
        String rawbam = read_string(stdout())
    }
}
# Sort reads by position
task bamsort {
    input {
        String bam
        String name = basename(bam)
        String dir
        String? temp = ""
        String temp_opt = if temp == "" then "" else "-T ~{temp}"
        Int? thread = 2
    }
    command <<<
        samtools sort -l1 \
          ~{temp_opt} \
          -@ ~{thread} \
          -O bam \
          -o ~{dir}/~{name}.sorted.bam \
          ~{bam}
        echo ~{dir}/~{name}.sorted.bam
    >>>
    output {
        String sorted = read_string(stdout())
    }
}
# Make BAM index (.bai)
task makeindex {
    input {
        String bam
    }
    command <<<
        samtools index ~{bam}
    >>>
    output {}
}