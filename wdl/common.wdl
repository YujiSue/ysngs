version 1.0
task curldl {
    input {
        String url
        String out
    }
    command <<< 
        curl -L -o ~{out} ~{url}
    >>>
    output {
        String dl = "~{out}"
    }
}
task getpath {
    input {
        String dir
        String out
        String ext
    }
    command <<<
        echo ~{dir}/~{out}.~{ext}
    >>>
    output {
        String path = read_string(stdout())
    }
}
task moveto {
    input {
        Array[String] ori
        String dest
    }
    command <<< 
        mv ~{sep=' ' ori} ~{dest}
        echo ~{dest}
    >>>
    output {
        String moved = read_string(stdout())
    }
}
task copyto {
    input {
        String from
        String to
    }
    command <<< 
        cp ~{from} ~{to}
    >>>
    output {}
}
task remove {
    input {
        Array[String] src
    }
    command <<<
        rm -r ~{sep=' ' src}
    >>>
    output {}
}

