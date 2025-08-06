nextflow.enable.dsl=2

process sayHi {
    output:
    file 'hello.txt'

    script:
    '''
    echo "Hi from Nextflow on SLURM!" > hello.txt
    '''
}

workflow {
    sayHi()
}
