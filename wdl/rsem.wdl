version 1.0
import "common.wdl" as common
task rsemcount1 {
	input {
		Array[File] fq
		String outdir
		String out
		String ref
		Int thread
	}
	command <<< 
		rsem-calculate-expression \
			-p ~{thread} \
			--hisat2-hca \
			--hisat2-path /content/MyApp/hisat2 \
			~{sep=" " fq}
			$HYM_REF/RSEM/~{ref} \
			$HYM_DATA/~{outdir}/~{out}
	>>>
	output {
	
	}
}
task rsemcount2 {
	input {
		Array[File] fq
		String outdir
		String out
		String ref
		Int thread
	}
	command <<< 
		rsem-calculate-expression \
			--paired-end \
			-p ~{thread} \
			--hisat2-hca \
			--hisat2-path /content/MyApp/hisat2 \
			~{sep=" " fq}
			$HYM_REF/RSEM/~{ref} \
			$HYM_DATA/~{outdir}/~{out}
	>>>
	output {
	
	}
}

task rsemindex {
	input {
		File fa
		File gtf
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/RSEM
		rsem-prepare-reference \
			-p ~{thread} \
			--gtf ~{gtf} \
			~{fa} \
			$HYM_REF/RSEM/~{label}
		echo $HYM_REF/RSEM/~{label}
	>>>
	output {
		String refdir = read_string(stdout())
	}
}
task rsemindexb {
	input {
		File fa
		File gtf
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/RSEM
		rsem-prepare-reference \
			--gtf ~{gtf} \
			--num-threads ~{thread} \
			~{fa} \
			$HYM_REF/RSEM/~{label}
		echo $HYM_REF/RSEM/~{label}
	>>>
	output {
		String refdir = read_string(stdout())
	}
}
task rsemindexs {
	input {
		File fa
		File gtf
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/RSEM
		rsem-prepare-reference \
			--gtf ~{gtf} \
			--num-threads ~{thread} \
			~{fa} \
			$HYM_REF/RSEM/~{label}
		echo $HYM_REF/RSEM/~{label}
	>>>
	output {
		String refdir = read_string(stdout())
	}
}
task rsemindexh {
	input {
		File fa
		File gtf
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/RSEM
		rsem-prepare-reference \
			-p ~{thread} \
			--hisat2-hca \
			--hisat2-path $HYM_APP/hisat2 \
			--gtf ~{gtf} \
			~{fa} \
			$HYM_REF/RSEM/~{label}
		echo $HYM_REF/RSEM/~{label}
	>>>
	output {
		String refdir = read_string(stdout())
	}
}