rule align_consensus:
    input:
        cons=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.fasta',
        ref=config['reference']['fasta'],
    output:
        bam=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.{ref}.bam',
        idx=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.{ref}.bam.bai',
    threads: 
        2
    resources:
        mem_mb=16000,
        disk_mb=16000,
    log:
        f'batches/{batch}/logs/pbmm2/align_consensus.{{sample}}.log'
    conda:
        'envs/pbmm2.yaml'
    shell:
        '''
        (pbmm2 align -j {threads} \
                     --preset hifi \
                     --sort \
                     {input.ref} {input.cons} {output.bam}) > {log} 2>&1
        '''

rule align_hifi:
    input:
        fq=f'batches/{batch}/{{sample}}/primertrim_subsample/{{sample}}.fastq',
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
