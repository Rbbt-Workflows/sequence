Functionalities regarding genomic sequence analysis.

Positions, ranges, and mutations are specified using a common format. Each is
identified by several fields separated by the character `:`. *Positions* are
represented as `chromosome:position` i.e. `12:4234412`. *Mutations* have an
additional field representing the mutant allele `chromosome:position:allele`
i.e. `12:4234412:T` (the reference allele is redundant, as it is specified by
the chromosome and position). Indels are represented as in the following
examples: `+A` or `+ATC` for one and three base insertions repectively or `-`
and `---` for one or three deletions repectively. *Chromosome ranges* are
specified as `chromosome:start:end` as in `12:4234412:4244412`.

It supports multiple organisms. The format of the *organism* input is the
organism short code (`Hsa` for `H`omo `sa`piens, or `Mmu` for `M`us `mu`sculus)
optionally followed by the date of the build. For example, `Hsa/jan2013` for a
recent build or `Hsa/may2009` for the hg18 build.

The *watson* input is used to specify if the variants are described in
reference to the watson or forward strand, or in reference to the strand that
holds the overlapping gene. Using the wrong convention may make some mutations
coincide with the reference. The *is_watson* method can take a guess by
checking this criteria.

