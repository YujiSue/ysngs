version 1.0
workflow dlfq {
    input {
        String data_id
        String out_dir
        Int thread
    }
    call sradump { 
        input: 
          sraid = data_id,
          dir = out_dir,
          thread = thread
    }  
}
task sradump {
    input {
        String sraid
        String dir
        Int thread
    }
    command <<< 
        mkdir -p ~{dir}
        $HYM_APP/sra/fasterq-dump ~{sraid} -O ~{dir} -e ~{thread}
    >>>
    output {}
}
