localrules: extract_primer_data, extract_cluster_barcodes

def _get_primer_reports( wildcards ):
    return [ f'batches/{batch}/{wildcards.sample}/primertrim/{barcode}/primerClipped.lima.report'
             for barcode in sample2bc[ wildcards.sample ] ]                   

def _get_lima_report( wildcards ):
    demuxdir = checkpoints.demux_ubam.get( **wildcards ).output[0]
    return f'{demuxdir}/demultiplex.lima.report'

rule extract_primer_data:
    input:
        clustered_readnames=f'batches/{batch}/{{sample}}/clustered_holes.txt',
        primer_reports=_get_primer_reports,
    output:
        temp(f'batches/{batch}/{{sample}}/clusters.primer.report'),
    threads:
        1
    log:
        f'batches/{batch}/logs/extract/{{sample}}.primer_report.log'
    shell:
        '''
        (grep -hf <(cut -d/ -f1,2 {input.clustered_readnames}) {input.primer_reports} > {output}) > {log} 2>&1
        '''

rule extract_cluster_barcodes:
    input:
        clustered_readnames=f'batches/{batch}/{{sample}}/clustered_holes.txt',
        lima_report=_get_lima_report,
    output:
        temp(f'batches/{batch}/{{sample}}/clusters.barcode.report'),
    threads:
        1
    shell:
        '''
        grep -f <(cut -d/ -f1,2 {input.clustered_readnames}) {input.lima_report} > {output}
        '''

rule get_cluster_qc:
    input:
        pbaa_info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        demux_report=f'batches/{batch}/{{sample}}/clusters.barcode.report',
        ptrim_report=f'batches/{batch}/{{sample}}/clusters.primer.report',
    output:
        f'batches/{batch}/{{sample}}/clusterQC.report.tsv',
    params:
        minFreq=0.02,
    threads:
        1
    log:
        f'batches/{batch}/logs/python/{{sample}}.cluster_ends.log'
    run:
        import pandas as pd
        cluster_info = pd.read_csv( input.pbaa_info, sep='\s', engine='python', usecols=[0,1,9], names=['read','target','cluster'], index_col=0 )
        cluster_info.index = cluster_info.index.str.rsplit('/',n=1).str[0]
        end_calls = pd.merge(
                            *[ pd.read_csv( tbl, sep='\t', usecols=[0,7,35,36], names=['read','qual','fwd','rev'], index_col=0 )
                            for tbl in [ input.demux_report, input.ptrim_report ] ],
                            left_index=True, right_index=True,
                            suffixes=[ '_barcode','_primer' ] 
                            ).join( cluster_info )
        counts = end_calls.groupby( ['target','cluster', 'fwd_primer','rev_primer'] ).describe(percentiles=[0.9]).sort_index()
        freqs = counts[ ('qual_barcode','count') ] / len(end_calls)
        freqs.name = ('','frequency')
        counts = counts.join( freqs ).set_index( ('','frequency'), append=True )
        counts.index.names = counts.index.names[:-1] + ['frequency']
        counts.query( 'frequency >= @params.minFreq' )\
              .to_csv( output[0], sep='\t', float_format='%.3f' )


extra_targets.append(
    lambda wildcards:
        [
         f'batches/{batch}/{sample}/clusterQC.report.tsv'
         for sample in { bc2sample[bc] for bc in _get_bam_demuxed( wildcards ) }
        ]
 )
