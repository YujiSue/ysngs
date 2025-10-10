version 1.0

# DL from GCS
task dlgcs {
    input {
        String url
        String out
    }
    command <<<
        gsutil cp ~{url} ~{out} 
    >>>
    output {
        String dl = out
    }
}

# Deepvariant variant call
task gdvcall {
    input {
        Boolean use_gpu = true
        Int gpu = 0
        String gopt = if use_gpu then "--gpus ~{gpu}" else ""
        String ver
        String path

        String reference
        String refdir = sub(reference, "/[^/]+$", "")
        String refname = basename(reference)

        Boolean target_seq = false
        String target
        String beddir = sub(target, "/[^/]+$", "")
        String bed = basename(target)
        String tset = if target_seq then "-v \"~{beddir}\":\"/TARGET_DIR\"" else ""
        String topt = if target_seq then "--regions /TARGET_DIR/~{bed}" else ""

        String bam
        String bamdir = sub(bam, "/[^/]+$", "")
        String in = basename(bam)

        String dir
        String name
        String outdir = "~{dir}/variant"
        String out1 = "/OUTPUT_DIR/~{dir}/variant/~{name}.g.vcf"
        String out2 = "/OUTPUT_DIR/~{dir}/variant/~{name}.vcf"
        
        Int thread
    }
    command <<<
        mkdir -p ~{outdir}
        docker run \
          ~{gopt} \
          -v "~{refdir}":"/REF_DIR" \
          -v "~{bamdir}":"/INPUT_DIR" \
          -v "~{outdir}":"/OUTPUT_DIR" \
          ~{tset} \
          google/deepvariant:~{ver} ~{path} \
          --model_type=WGS \
          --ref /REF_DIR/~{refname} \
          --reads /INPUT_DIR/~{in} \
          ~{topt} \
          --num_shards ~{thread} \
          --output_gvcf "~{out1}" \
          --output_vcf "~{out2}" 
    >>>
    output {
        String gvcf = out1
        String vcf = out2
    }
}