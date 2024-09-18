version 1.0

task hisatmap {
	input {
		Array[File] fq
		String outdir
		String out
		String ref
		Int thread
	}
	command <<<
		$HYM_APP/hisat2/hisat2 -x $HYM_REF/HISAT/~{ref} -U ~{fq[0]} -p ~{thread} -S $HYM_DATA/~{outdir}/~{out}.hisat.sam
		echo $HYM_DATA/~{outdir}/~{out}.hisat.sam
	>>>
	output {
		String sam = read_string(stdout())
	}
}
task hisatmap2 {
	input {
		Array[File] fq
		String outdir
		String out
		String ref
		Int thread
	}
	command <<< 
		$HYM_APP/hisat2/hisat2 -x $HYM_REF/HISAT/~{ref} -1 ~{fq[0]} -2 ~{fq[1]} -p ~{thread} -S $HYM_DATA/~{outdir}/~{out}.hisat.sam 
		echo $HYM_DATA/~{outdir}/~{out}.hisat.sam
	>>>
	output {
		String sam = read_string(stdout())
	}
}

task hisatexon {
	input {
		File gtf
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/HISAT
		mkdir -p $HYM_REF/HISAT/~{label}
		python $HYM_APP/hisat2/scripts/hisat2_extract_exons.py ~{gtf}
		
	>>>
	output {
		String exons = read_string(stdout())
	}
}
task hisatss {
	input {
		File gtf
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/HISAT
		mkdir -p $HYM_REF/HISAT/~{label}
		python $HYM_APP/hisat2/scripts/hisat2_extract_splice_sites.py ~{gtf}
		
	>>>
	output {
		String splices = read_string(stdout())
	}
}
task hisatindex {
	input {
		File fa
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/HISAT
		mkdir -p $HYM_REF/HISAT/~{label}
		$HYM_APP/hisat2/hisat2-build -p ~{thread} ~{fa} $HYM_REF/HISAT/~{label}
		echo $HYM_REF/HISAT/~{label}
	>>>
	output {
		String refdir = read_string(stdout())
	}
}
task hisatindext {
	input {
		File fa
		File gtf
		String outdir
		String label
		Int thread
	}
	command <<< 
		mkdir -p $HYM_REF/HISAT
		mkdir -p $HYM_REF/HISAT/~{outdir}
		python $HYM_APP/hisat2/scripts/hisat2_extract_exons.py ~{gtf}
		python $HYM_APP/hisat2/scripts/hisat2_extract_splice_sites.py ~{gtf}
		$HYM_APP/hisat2/hisat2-build -p ~{thread} ~{fa} $HYM_REF/HISAT/~{label}
		
		echo $HYM_REF/HISAT/~{label}
	>>>
	output {
		String refdir = read_string(stdout())
	}
}
