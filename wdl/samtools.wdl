version 1.0

# BAM => single-read FASTQ
task bam2fq {
    input {
        String dir
        Boolean paired
        String bam
        String name = basename(bam)
        String out = "~{dir}/~{name}.fq"
    }
    command <<<
        mkdir -p ~{dir}
        samtools fastq ~{bam} > ~{dir}/~{name}.fq
    >>>
    output {
        String fq = out
    }
}
# BAM => paired-end FASTQ
task bam2pair {
    input {
        String dir
        String bam
        String name = basename(bam)
        String out1 = "~{dir}/~{name}_1.fq"
        String out2 = "~{dir}/~{name}_2.fq"
    }
    command <<<
        mkdir -p ~{dir}
        samtools collate -u -O ~{bam} | samtools fastq -1 ~{out1} -2 ~{out2} -0 /dev/null -s /dev/null -n
    >>>
    output {
        Array[String] fq = [out1, out2]
    }
}
# Make fasta index
task makefaidx {
    input {
        String fa
    }
    command <<<
        samtools faidx ~{fa}
    >>>
    output {
        String fai = "~{fa}.fai"
    }
}
# SAM => BAM
task sam2bam {
    input {
        String sam
        String name = basename(sam)
        String dir
        String out = "~{dir}/~{name}.bam"
        Int thread = 2
    }
    command <<<
        samtools view -@ ~{thread} \
          -b -o ~{out} \
          ~{sam} 
    >>>
    output {
        String rawbam = out
    }
}
# Sort reads by position
task bamsort {
    input {
        String bam
        String dir
        String name = basename(bam)
        String out = "~{dir}/~{name}.bam"
        String temp = ""
        String temp_opt = if temp == "" then "" else "-T ~{temp}"
        Int thread = 2
    }
    command <<<
        samtools sort -l1 \
          ~{temp_opt} \
          -@ ~{thread} \
          -O bam \
          -o ~{out} \
          ~{bam}
    >>>
    output {
        String sorted = out
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
    output {
        String index = "~{bam}.bai"
    }
}
# BCFtools variant call
task bcfcall {
    input {
        String reference
        String bam
        String dir
        String name = basename(bam)
        String out = "~{dir}/~{name}.vcf"
    }
    command <<<
        mkdir -p ~{dir}
        bcftools mpileup -Ou \
          -f ~{reference} \
          ~{bam} | bcftools call -vm -Oz \
          -o ~{out}
    >>>
    output {
        String vcf = out
    }
}