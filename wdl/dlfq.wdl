version 1.0
# Download fastq from archive server
workflow dlfq {
    input {
        Boolean from_sra = true
        Boolean split_file = false
        String data_id
        String out_dir
        Int thread = 2
    }
     # From NCBI SRA
    if (from_sra) {
        call sradump { 
            input: 
                split = split_file,
                sraid = data_id,
                dir = out_dir,
                thread = thread
        }  
    }
    output {}
}
# DL via SRAToolkit
task sradump {
    input {
        Boolean split = false
        String opt = if split then "--split-files" else ""
        String sraid
        String dir
        Int thread
    }
    command <<< 
        mkdir -p ~{dir}
        $HYM_APP/sra/fasterq-dump ~{opt} ~{sraid} -O ~{dir} -e ~{thread}
    >>>
    output {}
}
