version 1.0
# single count
task sccount {
    input {
	String reference
        String bam
        String dir
	String subdir
	String in_dir = "~{dir}/~{subdir}"
        String prefix
	String outname
        String out_dir = "~{dir}/~{outname}"

        Int core = 0
        String core_opt = if core > 0 then "--localcores=~{core}" else ""
        Int ram = 0
        String ram_opt = if ram > 0 then "--localmem=~{ram}" else ""
    }
    command <<<
        $HYM_APP/cellranger/bin/cellranger count \
          --id=~{outname} \
          --transcriptome=~{reference} \
          --fastqs=~{in_dir} \
          --sample=~{prefix} \
          --create-bam=~{bam} \
          ~{core_opt} \
          ~{ram_opt}
        mv ~{outname} ~{out_dir}
    >>>
    output {
        String result = out_dir
    }
}
# multi count
task scmulti {
    input {
        String dir
        String config
        String outname
        String out_dir = "~{dir}/~{outname}"

        Int core = 0
        String core_opt = if core > 0 then "--localcores=~{core}" else ""
        Int ram = 0
        String ram_opt = if ram > 0 then "--localmem=~{ram}" else ""
    }
    command <<<
        $HYM_APP/cellranger/bin/cellranger multi \
          --id=~{outname} \
          --csv=~{config} \
          ~{core_opt} \
          ~{ram_opt}
        mv ~{outname} ~{out_dir}
    >>>
    output {
        String result = out_dir
    }
}
# make reference
task scmkref {
    input {
        String fasta
        String gtf
        String out
    }
    command <<<
        cellranger mkref \
          --genome=~{out} \
          --fasta=~{fasta} \
          --genes=~{gtf}
    >>>
    output {
        String reference = out
    }
}
