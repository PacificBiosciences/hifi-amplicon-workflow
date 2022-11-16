rule primer_trim:
    input: 
        primers=config['primers'],
        reads=f'batches/{batch}/demux/demultiplex.{{barcode}}.bam',
    output:
        f'batches/{batch}/{{sample}}/primertrim/{{barcode}}/primerClipped.bam',
    params:
        preset=f'--ccs --min-score {config["minPrimerMatch"]} --min-end-score {config["minPrimerMatch"]} --min-ref-span 0.75 --different --min-scoring-regions 2',
        loglevel='INFO',
    threads:
        2
    log:
        f'batches/{batch}/logs/primertrim/{{sample}}.{{barcode}}.log'
    conda:
        'envs/lima.yaml'
    shell:
        '''
        (lima {params.preset} \
              -j {threads} \
              --log-level {params.loglevel} \
              {input.reads} {input.primers} {output}) > {log} 2>&1
        '''

#rule bam_to_fastq:
#    input:
#        [ f'batches/{batch}/{{sample}}/primertrim/{barcode}/primerClipped.bam'
#    output:
#        f'batches/{batch}/{{sample}}/primertrim/primerClipped.fastq',
#    threads:
#        2
#    log:
#        f'batches/{batch}/logs/samtools/bam_to_fastq.{{sample}}.log'
#    conda:
#        'envs/samtools.yaml'
#    shell:
#        '''
#        (samtools fastq {input} > {output}) > {log} 2>&1
#        '''
#
#rule fqidx:
#    input:
#        f'batches/{batch}/{{sample}}/primertrim/primerClipped.fastq',
#    output:
#        f'batches/{batch}/{{sample}}/primertrim/primerClipped.fastq.fai',
#    threads:
#        1
#    log:
#        f'batches/{batch}/logs/samtools/fqidx.{{sample}}.log'
#    conda:
#        'envs/samtools.yaml'
#    shell:
#        '''
#        (samtools fqidx {input}) > {log} 2>&1
#        '''
#
