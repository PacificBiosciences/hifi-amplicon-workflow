rule extract_clustered_reads:
    input:
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        fq=f'batches/{batch}/{{sample}}/primertrim_subsample/{{sample}}.fastq'
        idx=f'batches/{batch}/{{sample}}/primertrim_subsample/{{sample}}.fastq.fai',
    output:
        incl=temp(f'batches/{batch}/{{sample}}/clustered_holes.txt'),
        fq=temp(f'batches/{batch}/{{sample}}/clustered_hifi.fastq'),
    threads:
        1
    conda:
        'envs/samtools.yaml'
    log:
        f'batches/{batch}/logs/samtools/extract.{{sample}}.log'
    shell:
        '''
        (cut -d' ' -f1 {input.info} > {output.incl}
         seqtk subsample {input.fq} {output.incl} > {output.fq}) > {log} 2>&1
        '''

rule align_hifi:
    input:
        fq=f'batches/{batch}/{{sample}}/clustered_hifi.fastq',
        ref=config['reference']['fasta'],
    output:
        bam=temp(f'batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.aligned.bam'),
        idx=temp(f'batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.aligned.bam.bai'),
    threads: 
        8
    log:
        f'batches/{batch}/logs/pbmm2/align_hifi.{{sample}}.log'
    conda:
        'envs/pbmm2.yaml'
    shell:
        '''
        (pbmm2 align -j {threads} \
                     --preset hifi \
                     --sort \
                     {input.ref} {input.fq} {output.bam}) > {log} 2>&1
        '''

rule paint_bam:
    input:
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        bam=f'batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.aligned.bam',
        idx=f'batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.aligned.bam.bai',
    output:
        f'batches/{batch}/{{sample}}/hifi.painted.bam'
    threads: 
        1
    conda:
        'envs/pbaa.yaml'
    log:
        f'batches/{batch}/logs/pbaa/bampaint.{{sample}}.log'
    shell:
        '''
        (pbaa bampaint {input.info} {input.bam} {output}) > {log} 2>&1
        '''

extra_targets.append( 
    lambda wildcards:
        [
         f'batches/{batch}/{sample}/hifi.painted.bam'
         for sample in { bc2sample[bc] for bc in _get_bam_demuxed( wildcards ) } 
        ]
 ) 
