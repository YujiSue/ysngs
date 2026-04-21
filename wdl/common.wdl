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
task makelist {
    input {
        String target
        String out
    }
    command <<<
        ls ~{target} > ~{out}
    >>>
    output {
        String list = out
    }
}
task concat {
    input {
        Array[String] files
        String out
    }
    command <<<
        cat ~{sep=" " files} > ~{out}
    >>>
    output {
        String concat = out
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

