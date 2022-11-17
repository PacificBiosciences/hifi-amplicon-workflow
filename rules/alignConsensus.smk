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
