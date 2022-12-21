rule primer_trim:
    input: 
        primers=config['primers'],
        reads=f'batches/{batch}/demux/demultiplex.{{barcode}}.bam',
    output:
        base=f'batches/{batch}/{{sample}}/primertrim/{{barcode}}/primerClipped.bam',
        report=f'batches/{batch}/{{sample}}/primertrim/{{barcode}}/primerClipped.lima.report',
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
              {input.reads} {input.primers} {output.base}) > {log} 2>&1
        '''
