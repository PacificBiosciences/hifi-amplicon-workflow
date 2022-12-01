rule report:
    input: 
        vcf=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.annotated.vcf.gz',
        vcf_index=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.annotated.vcf.gz.tbi',
    output:
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.variants.tsv',
    params:
        fields='\'%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%CLNSIG\t%CLNDN]\n\'',
    threads:
        1
    log:
        f"batches/{batch}/logs/bcftools/report_{{sample}}.{{consensus}}.log"
    conda:
        "envs/bcftools.yaml"
    message:
        "Reporting variants for {wildcards.sample}: {wildcards.consensus}"
    shell:
        '''
        (bcftools query -H -f {params.fields} {input.vcf} > {output}) > {log} 2>&1
        '''

def _get_alleles( wildcards ):
    extract_dir = checkpoints.extract_alignments.get( sample=wildcards.sample ).output[0]
    consensus   = glob_wildcards( f'{extract_dir}/{{consensus}}.{ref}.bam' ).consensus
    return expand( f'batches/{batch}/{wildcards.sample}/vcf/{{consensus}}.{ref}.htsbox.variants.tsv', consensus=consensus )

rule collate_alleles:
    input:
        tsvs=_get_alleles,
    output:
        f'batches/{batch}/{{sample}}/{{sample}}.{ref}.clinvar_annotated.variant_summary.tsv',
    run:
        import pandas as pd
        patt = 'sample-(?P<barcode>.*)_guide-(?P<guide>.*)_cluster-(?P<cluster>[0-9]+)_ReadCount-(?P<numreads>[0-9]+).*:(?P<field>.*$)'
        fmt = lambda r: f'{r.guide}_{r.cluster}_numreads{r.numreads}:{r.field}'
        
        def read_tsv( tsv ):
            tbl = pd.read_csv( tsv, sep='\t' )
            tbl.columns = tbl.columns.str.replace( '^.*\]', '', regex=True )
            tbl = tbl[tbl.ID != "."]
            infoCols = tbl.columns.str.contains( 'cluster' )
            tbl.columns = tbl.columns.where( ~infoCols, 
                                              tbl.columns.str.extract( patt ).apply( fmt, axis=1 )
                                            )
            return tbl.set_index( list( tbl.columns[ ~infoCols ] ) )
        pd.concat( map( read_tsv, sorted( input.tsvs ) ), axis=1 )\
          .dropna( axis=1, how='all' )\
          .fillna( '<not called>' )\
          .to_csv( output[0], sep='\t' )
