rule pbaa_make_guide:
    input:
        reference=config['reference']['fasta'],
    output:
        f'batches/{batch}/pbaa/guide.fasta',
        f'batches/{batch}/pbaa/guide.fasta.fai',
    params:
        regions=config['regions'],
    log:
        f'batches/{batch}/logs/pbaa/make_guide.log'
    run:
        import pysam
        with pysam.FastaFile( input.reference ) as infa, \
             open( output[0], 'w' ) as outfa:
            for label,region in params.regions.items():
                sequence = infa.fetch( region=region )
                outfa.write( f'>{label}\n{sequence}\n' )
        pysam.faidx( output[0] )


rule pbaa_cluster:
    input:
        fq=f'batches/{batch}/{{sample}}/primertrim_subsample/{{sample}}.fastq',
        idx=f'batches/{batch}/{{sample}}/primertrim_subsample/{{sample}}.fastq.fai',
        guide=f'batches/{batch}/pbaa/guide.fasta'
    output:
        cons1=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.fasta',
        cons2=f'batches/{batch}/{{sample}}/pbaa_failed_cluster_sequences.fasta',
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        log=f'batches/{batch}/{{sample}}/pbaa.log'
    params:
        prefix=f'batches/{batch}/{{sample}}/pbaa',
        loglevel='INFO',
        maxReads=config[ 'maxClusteringReads' ],
        minFreq=config[ 'minClusterFrequency' ],
        minReads=config[ 'minClusterReads' ],
        maxUchime=config[ 'maxUchime' ],
        maxLen=config[ 'maxAmpliconSize' ],
    threads:
        8
    log:
        f'batches/{batch}/logs/pbaa/cluster.{{sample}}.log'
    conda:
        'envs/pbaa.yaml'
    shell:
        '''
        (pbaa cluster --min-cluster-frequency {params.minFreq} \
                      --max-uchime-score {params.maxUchime} \
                      --max-reads-per-guide {params.maxReads} \
                      --min-cluster-read-count {params.minReads} \
                      --max-amplicon-size {params.maxLen} \
                      -j {threads} \
                      --log-file {output.log} \
                      --log-level {params.loglevel} \
                      {input.guide} {input.fq} {params.prefix}) > {log} 2>&1
        '''

#rule extract_clustered_reads:
#    input:
#        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
#        bam=f'batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.aligned.bam',
#        idx=f'batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.aligned.bam.bai',
#    output:
#        incl=temp(f'batches/{batch}/{{sample}}/clustered_holes.txt'),
#        tmp=temp(f'batches/{batch}/{{sample}}/hifi.bam'),
#        idx=temp(f'batches/{batch}/{{sample}}/hifi.bam.bai'),
#    threads:
#        1
#    conda:
#        'envs/samtools.yaml'
#    log:
#        f'batches/{batch}/logs/samtoole/extract.{{sample}}.log'
#    shell:
#        '''
#        (cut -d' ' -f1 {input.info} > {output.incl}
#         samtools view -bhN {output.incl} {input.bam} --write-index -o {output.tmp}##idx##{output.idx}) > {log} 2>&1
#        '''
#        
#rule paint_bam:
#    input:
#        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
#        bam=f'batches/{batch}/{{sample}}/hifi.bam',
#        idx=f'batches/{batch}/{{sample}}/hifi.bam.bai',
#    output:
#        f'batches/{batch}/{{sample}}/hifi.painted.bam'
#    threads: 
#        1
#    conda:
#        'envs/pbaa.yaml'
#    log:
#        f'batches/{batch}/logs/pbaa/bampaint.{{sample}}.log'
#    shell:
#        '''
#        (pbaa bampaint {input.info} {input.bam} {output}) > {log} 2>&1
#        '''

#rule index_bam:
#    input:
#        f'batches/{batch}/{{sample}}/hifi.painted.bam'
#    output:
#        f'batches/{batch}/{{sample}}/hifi.painted.bam.bai'
#    threads: 
#        1
#    log:
#        f'batches/{batch}/logs/samtools/index_painted.{{sample}}.log'
#    conda:
#        'envs/samtools.yaml'
#    shell: 
#        '''
#        (samtools index {input}) > {log} 2>&1
#        '''
