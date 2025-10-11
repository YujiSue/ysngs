version 1.0
import "common.wdl" as common
#
task rsemcount {
	input {
		Boolean include_mapping = true
		Boolean use_bowtie = false
		Boolean use_star = false
		Boolean use_hisat = false

		
		Boolean mapping = false
		String mapper = ""
		String mapper_path = ""
		String mapopt = ""

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
			~{type} \
			~{if mapping then mapopt} \
			~{if mapping then mapper_path} \
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
		Boolean include_mapping = true
		Boolean use_bowtie = false
		Boolean use_star = false
		Boolean use_hisat = false

		String opt = if include_mapping then (if use_bowtie then "--bowtie2" else (if use_star then "--star" else (if use_hisat then "--hisat2-hca" else ""))) else ""
		String mapper = if include_mapping then (if use_bowtie then "--bowtie2-path $HYM_APP/bowtie2" else (if use_star then "--star-path $HYM_APP/STAR/bin/Linux_x86_64" else (if use_hisat then "--hisat2-path $HYM_APP/hisat2" else ""))) else ""

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