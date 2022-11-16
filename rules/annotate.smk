rule annotate:
    input: 
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf',
    output:
        f'batches/{batch}/{{sample}}/vcf/{{consensus}}.{ref}.htsbox.vcf'
