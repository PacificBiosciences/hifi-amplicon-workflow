localrules: tabix_vcf, extract_alignments

checkpoint extract_alignments:
    input:
        bam=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.{ref}.bam',
        idx=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.{ref}.bam.bai',
    output:
        directory( f'batches/{batch}/{{sample}}/vcf' ),
    threads:
        1
    run:
        import os, pysam
        os.makedirs( output[0] )
        with pysam.AlignmentFile( input.bam ) as inbam:
            for rec in inbam:
                obam = f'{output[0]}/{rec.query_name}.{ref}.bam'
                pysam.AlignmentFile( obam, 'wb', header=inbam.header ).write( rec )
                pysam.index( obam )

rule htsbox:
    input:
        bam=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.bam',
        idx=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.bam.bai',
        reference = config['reference']['fasta'],
    output: 
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf'
    log: 
        f"batches/{batch}/logs/htsbox/{{sample}}.{{consensus}}.log"
    conda: 
        "envs/htsbox.yaml"
    message: 
        "Calling variants from {input.bam} using htsbox."
    shell: 
        '''
        (htsbox pileup -c -f {input.reference} {input.bam} > {output})> {log} 2>&1
        '''

def _agg_extracted_alignments( wildcards ):
    extract_dir = checkpoints.extract_alignments.get( **wildcards ).output[0]
    consensus   = glob_wildcards( f'{extract_dir}/{{consensus}}.{ref}.bam' ).consensus
    return expand( f'{extract_dir}/{{consensus}}.{ref}.bam',
                   sample=wildcards.sample, consensus=consensus )

rule bgzip_vcf:
    input:
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf'
    output:
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf.gz'
    log:
        f"batches/{batch}/logs/bgzip/{{sample}}.{{consensus}}.{ref}.htsbox.log",
    threads: 2
    conda:
        "envs/htslib.yaml"
    message:
        "Executing {rule}: Compressing {input}."
    shell:
        "(bgzip --threads {threads} {input}) > {log} 2>&1"

rule tabix_vcf:
    input:
        f'batches/{batch}/{{sample}}/{{prefix}}.vcf.gz'
    output:
        f'batches/{batch}/{{sample}}/{{prefix}}.vcf.gz.tbi'
    threads: 1
    conda:
        "envs/htslib.yaml"
    message:
        "Executing {rule}: Indexing {input}."
    shell:
        "tabix -p vcf {input}"
