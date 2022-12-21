rule annotate:
    input: 
        vcf=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf.gz',
        vcf_index=f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf.gz.tbi',
        reference=config["reference"]["fasta"],
        clinvar=config["annotation"]["clinvar"],
        conseq=config["annotation"]["consequence"],
    output:
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.annotated.vcf.gz',
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
                  -a {input.clinvar} \
                  {input.vcf} | \
         bcftools csq \
                  -f {input.reference} \
                  -g {input.conseq} | \
         bcftools reheader \
                  -s <(echo {wildcards.consensus}) | \
         bcftools view \
                  -Oz \
                  -o {output}) > {log} 2>&1
        '''
