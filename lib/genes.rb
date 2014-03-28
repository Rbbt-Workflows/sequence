module Sequence

  desc "Find genes at particular positions in a chromosome. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :positions, :array, "Positions"
  def self.genes_at_chr_positions(organism, chromosome, positions)
    index = gene_position_index(organism, chromosome)
    index.chunked_values_at(positions).collect{|list| list * "|"}
  end
  task :genes_at_chr_positions => :array
  export_exec :genes_at_chr_positions

  desc "Find genes at particular genomic positions. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 11:533766). Separator can be ':', space or tab. Extra fields are ignored"
  def self.genes_at_genomic_positions(organism, positions)
    raise ParameterException, "No positions given" if positions.nil?
    chr_positions = {}
    positions.each do |position|
      next if position.empty?
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    chr_genes = {}
    chr_positions.each do |chr, list|
      chr_genes[chr] = genes_at_chr_positions(organism, chr, list)
    end

    tsv = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism)
    positions.collect do |position|
      next if position.empty?
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      tsv[position] = chr_genes[chr].shift.split("|")
    end
    tsv
  end
  task :genes_at_genomic_positions => :tsv
  export_synchronous :genes_at_genomic_positions

  desc "Find genes close to particular genomic positions. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 11:533766). Separator can be ':', space or tab. Extra fields are ignored"
  input :upstream, :integer, "Upstream bases", 1000
  input :downstream, :integer, "Downstream bases", 100
  def self.genes_close_to_genomic_positions(organism, positions, upstream, downstream)
    margin = [upstream, downstream].max
    chr_positions = {}
    positions.each do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    chr_genes_forward = {}
    chr_genes_reverse = {}
    chr_positions.each do |chr, list|
      forward_ranges = list.collect{|position| [position.to_i - upstream, position.to_i + downstream] * ":" }
      reverse_ranges = list.collect{|position| [position.to_i - downstream, position.to_i + upstream] * ":" }
      chr_genes_forward[chr] = genes_at_chr_ranges(organism, chr, forward_ranges)
      chr_genes_reverse[chr] = genes_at_chr_ranges(organism, chr, reverse_ranges)
    end

    all_genes = chr_genes_forward.values.flatten.uniq + chr_genes_reverse.values.flatten.uniq
    all_genes.uniq!

    tsv = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism)
    strand_tsv = gene_strand_index(organism)
    positions.collect do |position|
      chr, pos = position.split(/[\s:\t]/).values_at 0, 1
      chr.sub!(/chr/,'')
      genes = chr_genes_forward[chr].shift.split("|").select{|gene| strand_tsv[gene] == "1"} + 
        chr_genes_reverse[chr].shift.split("|").select{|gene| strand_tsv[gene] == "-1"} 
      tsv[position] = genes
    end
    tsv
  end
  task :genes_close_to_genomic_positions => :tsv
  export_synchronous :genes_close_to_genomic_positions

  desc "Find genes at particular ranges in a chromosome. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :ranges, :array, "Ranges"
  def self.genes_at_chr_ranges(organism, chromosome, ranges)
    index = gene_position_index(organism, chromosome)
    r = ranges.collect{|r| s,e = r.split(":"); (s.to_i..e.to_i)}
    index.values_at(*r).collect{|list| list * "|"}
  end
  task :genes_at_chr_ranges => :array
  export_exec :genes_at_chr_ranges

  desc "Find genes at particular genomic ranges. Multiple values separated by '|'"
  input :organism, :string, "Organism code", "Hsa"
  input :ranges, :array, "Positions Chr:Start:End (e.g. 11:533766:553323)"
  def self.genes_at_genomic_ranges(organism, ranges)
    chr_ranges = {}
    ranges.each do |range|
      chr, s, e = range.split(":").values_at 0, 1, 2
      chr.sub!(/chr/,'')
      chr_ranges[chr] ||= []
      chr_ranges[chr] << [s, e] * ":" 
    end

    chr_genes = {}
    chr_ranges.each do |chr, list|
      chr_genes[chr] = genes_at_chr_ranges(organism, chr, list)
    end

    tsv = TSV.setup({}, :key_field => "Genomic Range", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism, :unnamed => true)
    ranges.collect do |range|
      chr, s, e = range.split(":").values_at 0, 1, 2
      chr.sub!(/chr/,'')
      tsv[range] = chr_genes[chr].shift.split("|")
    end
    tsv
  end
  task :genes_at_genomic_ranges => :tsv
  export_synchronous :genes_at_genomic_ranges
end
