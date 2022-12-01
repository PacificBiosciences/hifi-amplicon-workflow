rule annotate:
    input: 
        vcf=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf.gz',
        vcf_index=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf.gz.tbi',
        annot=config["annotation"],
    output:
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.annotated.vcf.gz'
    threads:
        1
    log:
        f"batches/{batch}/logs/bcftools/annotate_{{sample}}.{{consensus}}.log"
    conda:
        "envs/bcftools.yaml"
    message:
        "Annotating variants for {wildcards.sample}: {wildcards.consensus}"
    shell:
        '''
        (bcftools annotate \
                  -c ID,INFO \
                  -a {input.annot} \
                  -Oz \
                  -o {output} \
                  {input.vcf}) > {log} 2>&1
        '''
