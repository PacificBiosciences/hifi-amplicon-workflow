localrules: tabix_vcf

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

rule extract_clustered_reads:
    input:
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        fq=f'batches/{batch}/{{sample}}/primertrim_subsample/{{sample}}.fastq',
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
         seqtk subseq {input.fq} {output.incl} > {output.fq}) > {log} 2>&1
        '''
