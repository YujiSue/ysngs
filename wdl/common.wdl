version 1.0

task curldl {
    input {
        String url
        String out
    }
    command <<< 
        curl -L -o "~{out}" "~{url}"
    >>>
    output {
        String dl = out
    }
}
task copyto {
    input {
        Boolean copydir = false
        String opt = if copydir then "-r" else "" 
        String src
        String dest
    }
    command <<< 
        cp ~{opt} ~{src} ~{dest}
    >>>
    output {
        String to = dest
    }
}
task remove {
    input {
        Array[String] src
        Array[String] exclude
    }
    command <<<
        rm -r ~{sep=' ' src}
    >>>
    output {}
}

