from collections import defaultdict

shell.prefix(
     f"set -o pipefail; umask 002; export TMPDIR={config['tmpdir']};"
)

batch     = config[ 'batch' ]
hifiReads = config[ 'hifiReads' ]
ref       = config['reference']['label']

# sample-barcode map
# requires barcodeCsv to have two columns "barcode","sample" with a header line
# (header names not important, just the column order barcodes 1st, sample 2nd)
bc2sample = dict( 
                  line.strip().split( ',' ) 
                  for i,line in enumerate( open( config[ 'biosamples' ] ).readlines() )
                  if i > 0
                ) 
# Reverse mapping uses a list so more than one barcode can encode a single sample
sample2bc = defaultdict(list)
for k,v in bc2sample.items():
    sample2bc[v].append(k)
# add in barcodes found when samples have hybrids and deletions
if config[ 'demuxCross' ]:
    from itertools import product
    for sample, barcodes in sorted( sample2bc.items() ):
        for f,r in product( *zip( *[ b.split('--') for b in barcodes ] ) ):
            bc = f'{f}--{r}'
            sample2bc[sample].append(bc)
            bc2sample[bc] = sample

# write biosamples csv for use in lima
with open( f'batches/{batch}/biosamples.csv', 'w' ) as biosamples:
    biosamples.write( 'Barcode,BioSample\n' )
    for barcode,sample in bc2sample.items():
        biosamples.write( f'{barcode},{sample}\n' ) 

def _get_bam_demuxed( wildcards ):
    '''
    Some samples may not have yield (failed), so update samples after demuxing
    '''
    demuxdir = checkpoints.demux_ubam.get( **wildcards ).output[0]
    return glob_wildcards( f'{demuxdir}/demultiplex.{{barcode}}.bam').barcode

def _agg_consensus( wildcards ):
    '''
    Returns a list of expected consensus results, given samples that produced demux output
    '''
    samples = { bc2sample[bc] for bc in _get_bam_demuxed( wildcards ) }
    return expand( f'batches/{batch}/{{sample}}/pbaa_{{status}}_cluster_sequences.fasta',
                   sample=samples,
                   status=['passed','failed'] )

def _agg_mapped_consensus( wildcards ):
    '''
    Returns a list of expected mapped consensus results, given samples that produced demux output
    '''
    samples = { bc2sample[bc] for bc in _get_bam_demuxed( wildcards ) }
    return expand( f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.{ref}.bam', sample=samples )

def _agg_variant_summary( wildcards ):
    samples = { bc2sample[bc] for bc in _get_bam_demuxed( wildcards ) }
    return expand( f'batches/{batch}/{{sample}}/{{sample}}.{ref}.clinvar_annotated.variant_summary.xlsx', sample=samples )

extra_targets = []

include: "rules/common.smk"
include: "rules/demux.smk"
include: "rules/primertrim.smk"
include: "rules/bam2fastq.smk"
include: "rules/pbaa.smk"
include: "rules/alignConsensus.smk"
include: "rules/htsbox.smk"
include: "rules/annotate.smk"
include: "rules/report.smk"
include: "rules/clusterMetrics.smk"

if config['alignHiFi']:
    include: "rules/alignHiFi.smk"

rule all:
    input:
        _agg_consensus,
        _agg_mapped_consensus,
        _agg_variant_summary,
        extra_targets,
