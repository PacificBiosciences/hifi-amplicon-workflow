rule report:
    input: 
        vcf=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.annotated.vcf.gz',
        vcf_index=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.annotated.vcf.gz.tbi',
    output:
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.variants.tsv',
    params:
        fields=lambda wc: '\'%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%CLNSIG\t%CLNDN\t%TBCSQ{0}]\n\'',
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
        (bcftools query -Hu -f {params.fields} {input.vcf} > {output}) > {log} 2>&1
        '''

def _get_alleles( wildcards ):
    extract_dir = checkpoints.extract_alignments.get( sample=wildcards.sample ).output[0]
    consensus   = glob_wildcards( f'{extract_dir}/{{consensus}}.{ref}.bam' ).consensus
    return expand( f'batches/{batch}/{wildcards.sample}/vcf/{{consensus}}.{ref}.htsbox.variants.tsv', consensus=consensus )

rule collate_alleles:
    input:
        tsvs=_get_alleles,
    output:
        f'batches/{batch}/{{sample}}/{{sample}}.{ref}.clinvar_annotated.variant_summary.xlsx',
    run:
        import pandas as pd
        patt = 'sample-(?P<Sample>.*)_guide-(?P<Target>.*)_cluster-(?P<Cluster>[0-9]+)_ReadCount-(?P<Numreads>[0-9]+)'

        def read_tsv( tsv ):
            tbl = pd.read_csv( tsv, sep='\t' )
            tbl.columns = tbl.columns.str.replace( '^.*\]', '', regex=True ).str.split(':').str[-1]
            if tbl.empty:
                tbl.loc[0] = {'SAMPLE' : tsv, 'CHROM': 'NoVariants'}
            #fill in columns with missing cluster ("sample") name
            tbl.SAMPLE = tbl.SAMPLE[ tbl.SAMPLE.notnull()].iloc[0]
            clusterInfo = tbl.SAMPLE.str.extract( patt )
            res = pd.concat( [ tbl, clusterInfo ], axis=1 ).drop( columns='SAMPLE' ).fillna( '.' )
            res.set_index( [ 'Sample','Target','Cluster','Numreads','CHROM','POS'], inplace=True )
            return res
        
        pd.concat( map( read_tsv, input.tsvs ), axis=0 )\
          .sort_index()\
          .fillna('.')\
          .to_excel( output[0] )
