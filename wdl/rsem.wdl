version 1.0
import "common.wdl" as common
#
task rsemcount {
    input {
        Boolean use_bowtie = false
        Boolean use_star = false
        Boolean use_hisat = false
        Boolean export_bam = true
        String interopt = if export_bam then "--keep-intermediate-files --output-genome-bam" else ""
        
        Boolean mapping = false
        String opt = if mapping then (if use_bowtie then "--bowtie2" else (if use_star then "--star" else (if use_hisat then "--hisat2-hca" else ""))) else "--alignments"
        String mapper = if mapping then (if use_bowtie then "--bowtie2-path $HYM_APP/bowtie2" else (if use_star then "--star-path $HYM_APP/STAR/bin/Linux_x86_64" else (if use_hisat then "--hisat2-path $HYM_APP/hisat2" else ""))) else ""

        Array[String] fq = []
        String bam = ""
        Boolean paired
        String type = if paired then "--paired-end" else ""
        String dir
        String name
        String ref
        String out = "~{dir}/count/~{name}"

        Int thread = 2
    }
    command <<< 
        mkdir -p ~{dir}/count
        rsem-calculate-expression \
            -p ~{thread} \
            ~{interopt} \
            ~{type} \
            ~{opt} \
            ~{mapper} \
            ~{bam}~{sep=" " fq} \
            ~{ref} \
            ~{out}
        rm -r ~{out}.temp
    >>>
    output {
        String tbam = "~{out}.transcript.bam"
        String gbam = "~{out}.genome.bam"
        String count = "~{out}.genes.results"
        String count2 = "~{out}.isoforms.results"
    }
}
#
task rsemindex {
    input {
        Boolean mapping = true
        Boolean use_bowtie = false
        Boolean use_star = false
        Boolean use_hisat = false
		String mapper_path = ''

        String opt = if mapping then (if use_bowtie then "--bowtie2" else (if use_star then "--star" else (if use_hisat then "--hisat2-hca" else ""))) else ""
        String mapper = if mapping then (if use_bowtie then "--bowtie2-path ~{mapper_path}" else (if use_star then "--star-path ~{mapper_path}" else (if use_hisat then "--hisat2-path ~{mapper_path}" else ""))) else ""

        String fa
        String gtf
        String dir
        String label
        Int thread
    }
    command <<< 
        mkdir -p ~{dir}
        rsem-prepare-reference \
          --gtf ~{gtf} \
          --num-threads ~{thread} \
          ~{opt} \
          ~{mapper} \
          ~{fa} \
          ~{dir}/~{label}
    >>>
    output {
    }
}