version 1.0
import "common.wdl" as common
#
task rsemcount {
	input {
		Boolean use_bowtie = false
		Boolean use_star = false
		Boolean use_hisat = false
		Boolean export_bam = true
		String interopt = if export_bam then "--keep-intermediate-files" else ""
		
		Boolean mapping = false
		String opt = if mapping then (if use_bowtie then "--bowtie2" else (if use_star then "--star" else (if use_hisat then "--hisat2" else ""))) else ""
		String mapper = if mapping then (if use_bowtie then "--bowtie2-path $HYM_APP/bowtie2" else (if use_star then "--star-path $HYM_APP/STAR/bin/Linux_x86_64" else (if use_hisat then "--hisat2-path $HYM_APP/hisat2" else ""))) else ""

		Array[String] fq
		Boolean paired
		String type = if paired then "--paired-end" else ""
		String dir
		String name
		String ref

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
			~{sep=" " fq} \
			~{ref} \
			~{dir}/count/~{name}
	>>>
	output {
	
	}
}
#
task rsemindex {
	input {
		Boolean mapping = true
		Boolean use_bowtie = false
		Boolean use_star = false
		Boolean use_hisat = false

		String opt = if mapping then (if use_bowtie then "--bowtie2" else (if use_star then "--star" else (if use_hisat then "--hisat2-hca" else ""))) else ""
		String mapper = if mapping then (if use_bowtie then "--bowtie2-path $HYM_APP/bowtie2" else (if use_star then "--star-path $HYM_APP/STAR/bin/Linux_x86_64" else (if use_hisat then "--hisat2-path $HYM_APP/hisat2" else ""))) else ""

		String fa
		String gtf
		String dir
		String label
		Int thread
	}
	command <<< 
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