require 'rbbt/util/R'

module Sequence

  dep :affected_genes
  input *MUTATIONS_INPUT
  task :binomial_significance => :tsv do |mutations|
    mutation_count = Array === mutations ? Misc.counts(mutations.collect{|m| m.split(":").values_at(0,1)*":"}) : Misc.counts(Open.read(mutations).split("\n").collect{|m| m.split(":").values_at(0,1)*":"})
    iif mutation_count
    genes_step = step(:affected_genes)
    organism = genes_step.step(:mutated_isoforms_fast).inputs[1]

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Matches", "Bases", "Frequency"], :namespace => organism, :type => :list, :cast => :to_f)

    genes = Set.new
    mutations = Set.new
    gene_mutations = {}
    TSV.traverse genes_step do |position,g|
      position = position.split(":").values_at(0,1)*":"
      mutations << position
      next if g.nil?
      g.each do |gene|
        genes << gene
        gene_mutations[gene] ||= []
        gene_mutations[gene] << position
      end
    end

    genes = genes.to_a

    #num_mutations = mutations.to_a.length
    num_mutations = mutation_count.values.inject(0){|acc,e| acc += e; acc}

    total_bases = Organism.gene_list_exon_bases(genes)
    global_frequency = num_mutations.to_f / total_bases

    gene2exon_size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_exon_bases([gene]) }}

    TSV.traverse genes, :bar => "Computing significance" do |gene|
      mutations = gene_mutations[gene].uniq
      next if mutations.empty?
      matches = mutations.collect{|mutation| mutation_count[mutation]}.inject(0){|acc,e| acc += e}
      exon_bases = gene2exon_size[gene]
      next if exon_bases == 0
      frequency = matches.to_f / exon_bases

      tsv[gene] = [matches, exon_bases, frequency]
    end

    begin
      new = tsv.R "
      data = cbind(data,p.value=apply(data, 1, function(v){v = as.numeric(v); binom.test(v[1], v[2], #{ global_frequency }, 'greater')$p.value}))
      ", :key => "Ensembl Gene ID" 
      new.namespace = organism
      new
    rescue Exception
      tsv.annotate({}).add_field "p.value" do
      end
    end
  end
  export_asynchronous :binomial_significance

  dep :genes
  task :binomial_significance_syn => :tsv do 
    genes_step = step(:genes)
    organism = genes_step.info[:inputs][:organism]

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Matches", "Bases", "Frequency"], :namespace => organism, :type => :list, :cast => :to_f)

    genes = Set.new
    mutations = Set.new
    gene_mutations = {}
    TSV.traverse genes_step do |position,g|
      position = position.split(":").values_at(0,1)*":"
      mutations << position
      next if g.nil?
      g.each do |gene|
        genes << gene
        gene_mutations[gene] ||= []
        gene_mutations[gene] << position
      end
    end

    genes = genes.to_a
    num_mutations = mutations.to_a.length

    total_bases = Organism.gene_list_exon_bases(genes)
    global_frequency = num_mutations.to_f / total_bases

    gene2exon_size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_exon_bases([gene]) }}

    TSV.traverse genes, :bar => "Computing significance" do |gene|
      mutations = gene_mutations[gene].uniq
      next if mutations.empty?
      matches = mutations.length
      exon_bases = gene2exon_size[gene]
      next if exon_bases == 0
      frequency = matches.to_f / exon_bases

      tsv[gene] = [matches, exon_bases, frequency]
    end

    new = tsv.R "
      data = cbind(data,p.value=apply(data, 1, function(v){v = as.numeric(v); binom.test(v[1], v[2], #{ global_frequency }, 'greater')$p.value}))
      data
    ", :key => "Ensembl Gene ID" 

    new.namespace = organism
    new
  end
  export_asynchronous :binomial_significance_syn
end

