module Sequence

  input :positions, :array, "Positions Chr:Position (e.g. 19:54646887). Extra fields are ignored but kept", nil, :stream => false
  task :binomial_significance => :tsv do |positions|
    begin
      require 'rsruby'
    rescue
      raise "You must install rsruby gem to run this task"
    end

    genes_step = Sequence.job(:genes, name, :positions => positions)
    position_genes = genes_step.run
    organism = genes_step.info[:inputs][:organism]
    num_mutations = genes_step.info[:input_size] || positions.length

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Matches", "Bases", "Frequency", "p.value"], :namespace => organism, :type => :list, :cast => :to_f)

    genes = Set.new
    gene_mutations = {}
    TSV.traverse position_genes do |position,g|
      next if g.nil?
      g.each do |gene|
        genes << gene
        gene_mutations[gene] ||= []
        gene_mutations[gene] << position
      end
    end

    genes = genes.to_a

    total_bases = Organism.gene_list_exon_bases(genes)
    global_frequency = num_mutations.to_f / total_bases

    gene2exon_size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_exon_bases([gene]) }}

    genes.each do |gene|
      mutations = gene_mutations[gene]
      next if mutations.empty?
      matches = mutations.length
      exon_bases = gene2exon_size[gene]
      next if exon_bases == 0
      frequency = matches.to_f / exon_bases
      ddd  [matches, exon_bases, global_frequency, 'greater', gene]
      iii mutations
      pvalue = RSRuby.instance.binom_test(matches, exon_bases, global_frequency, 'greater')["p.value"]
      tsv[gene] = [matches, exon_bases, frequency, pvalue]
    end

    tsv
  end
  export_asynchronous :binomial_significance
end

