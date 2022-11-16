rule bam_to_fastq:
    input:
        lambda wcs:
            { f'batches/{batch}/{wcs.sample}/primertrim/{barcode}/primerClipped.bam'
             for barcode in sample2bc[wcs.sample] },
    output:
        f'batches/{batch}/{{sample}}/primertrim_subsample/{{sample}}.fastq',
    params:
        size=config['subsample'],
        seed=config[ 'randomSeed' ],
    threads:
        2
    log:
        f'batches/{batch}/logs/samtools/bam_to_fastq.{{sample}}.log'
    conda:
        'envs/samtools.yaml'
    shell:
        '''
        (parallel -j {threads} 'seqtk sample -s {params.seed} <(samtools fastq {{}}) {params.size}' ::: {input} > {output}) > {log} 2>&1
        '''

rule fqidx:
    input:
        f'batches/{batch}/{{prefix}}.fastq',
    output:
        f'batches/{batch}/{{prefix}}.fastq.fai',
    threads:
        1
    conda:
        'envs/samtools.yaml'
    shell:
        '''
        samtools fqidx {input}
        '''
