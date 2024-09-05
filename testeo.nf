log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}


process CONCATENATE{
    tag "CONCATENATION"
    input:
    tuple val(sample), val (order), val(lane), path(files)
    output:
    path "${sample}_${order}.fq"
    script:
    """
    echo "uniendo las partes del $order de $sample, $files" > "${sample}_${order}.fq"
    """
}

process TRIM{
    input:
    tuple val(sample), path(file)
    output:
    tuple val(sample), 
        path("${sample}_trimmed_R1.fq"),
        path("${sample}_trimmed_R2.fq")
    script:
    """
    echo "corto ${file[0]}" > "${sample}_trimmed_R1.fq"
    echo "corto ${file[1]}" > "${sample}_trimmed_R2.fq"
    """
}

process CUANTI {
    input:
    tuple val(sample), path(file1), path(file2)
    path salmon_index
    output:
    stdout
    script:
    """
     salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 $file -2 $file2 -o $sample
    """
}

workflow {
    Channel.fromPath("$projectDir/sampledata/*.fq")
    .map { file -> 
        def name = file.name
        def parts = name.tokenize('_')
        [parts[0], parts[2].replace('.fq', ''), parts[1], file]
    }
    .groupTuple(by:[0,1])
    | set {archivos_ch}
    CONCATENATE(archivos_ch)
    | map {it ->
    def name = it.name
    def parts = name.toString().tokenize('_')
    [parts[0],  it]
    }
    | groupTuple
    | map{sample, files ->
        sortedfiles = files.sort { file ->
            file.name.tokenize('_')[1].replace('.fq', '') == 'R2' ? 1 : 0
        }
        [sample,sortedfiles]
    }
    | TRIM
    | CUANTI


}