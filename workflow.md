Functionalities regarding genomic sequence analysis.

Finds genomic features overlapping genomic positions, like exons, reconstructs
offsets into transcripts, and computes the amino-acid changaes of variants.
Additionally finds mutations in exon junctions, and genes with high frequencies
of mutations.

# Tasks

Positions, ranges, and mutations are specified using a common format. Each is
identified by several fields separated by the character `:`. `Positions` are
represented as `chromosome:position` i.e. `12:4234412`. `Mutations` have an
additional field representing the mutant allele `chromosome:position:allele`
i.e. `12:4234412:T` (the reference allele is redundant, as it is specified by
the chromosome and position). Indels are represented as in the following
examples: `+A` or `+ATC` for one and three base insertions repectively or `-`
and `---` for one or three deletions repectively. `Chromosome ranges` are
specified as `chromosome:start:end` as in `12:4234412:4244412`.

It supports multiple organisms. The format of the `organism` input is the
organism short code (`Hsa` for `H`omo `sa`piens, or `Mmu` for `M`us `mu`sculus)
optionally followed by the date of the build. For example, `Hsa/jan2013` for a
recent build or `Hsa/may2009` for the hg18 build.

The `watson` input is used to specify if the variants are described in
reference to the watson or forward strand, or in reference to the strand that
holds the overlapping gene. Using the wrong convention may make some mutations
coincide with the reference. The `is_watson` method can take a guess by
checking this criteria.

Specifying the `vcf` parameter will interpret the input as a VCF file, and will
run the `genomic_mutations` task to extract the mutations from it

## reference

Report the reference base at the provided positions

## gene_strand_reference

Report the reference base at the provided positions on the gene coding strand

The gene coding strand is determined by checking for overlaping transcripts at
that position. In case of overlap the forward or watson strand is used.

## is_watson

Guess wether the mutations provided are given in the watson strand or the gene
strand

## genes

Report genes overlapping positions

## exons

Report exons overlapping positions

## transcripts

Report transcripts overlapping positions

## exon_junctions

Report exon junctions overlapping positions

## genes_at_ranges

Report genes overlapping ranges

## type

Report the type of base change: transition, transversion, indel, unknown or
nont at all

## transcript_offsets

Computes the offset inside the coding sequence of the transcripts overlapped the genomic mutations that
overlap them. 

## mutated_isoforms

Computes the consequence of genomic mutations in terms of amino-acid changes in
protein isoforms

All isoforms of a gene are reported

## splicing_mutations

Find mutations that may affect the splicing of protein coding transcripts

## mutated_isoforms_fast

One-step implementation of the `mutated_isoforms` task

## affected_genes

Finds genes affeted by genomic mutations, either by amino-acid changes on their
protein products, or by changes in splicing sequences

## binomial_significance

For a list of mutations, find genes that suffer a higher rate of mutation than
expected

Uses a binomial distribution with a global probablity of mutation estimated
from the data. Considers only exon bases

## expanded_vcf

Expands the `INFO` and `FORMAT/Sample` fields of VCF files in to a standard TSV format

## genomic_mutations

Extract genomic mutations from a VCF file that match a quality criteria

