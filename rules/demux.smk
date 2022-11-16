checkpoint demux_ubam:
    input:
        ccs=hifiReads,
        barcodes=config[ 'barcodes' ],
        biosamples=f'batches/{batch}/biosamples.csv'
    output:
        directory(f'batches/{batch}/demux')
    params:
        preset='--hifi-preset ASYMMETRIC --split-named',
        loglevel='INFO',
    threads:
        16
    log:
        f'batches/{batch}/logs/demux/lima.log'
    conda:
        'envs/lima.yaml'
    shell:
        '''
        (mkdir -p {output}
         lima {params.preset} \
              -j {threads} \
              --log-level {params.loglevel} \
              --biosample-csv {input.biosamples} \
              {input.ccs} {input.barcodes} {output}/demultiplex.bam) > {log} 2>&1
        '''

rule check_demux_fail:
    input:
        f'batches/{batch}/demux/',
    params:
        barcodes=_get_bam_demuxed,
    log:
        f"batches/{batch}/logs/lima/demux_no_yield.log",
    run:
        missing = bc2sample.keys() - params.barcodes
        if len( missing ):
            with open( f"batches/{batch}/demux_no_yield.txt", 'w' ) as ofile:
                for sample in missing:
                    ofile.write( f'{sample2barcode[sample]},{sample}\n' ) 

targets.append( f'batches/{batch}/logs/lima/demux_no_yield.log' )
