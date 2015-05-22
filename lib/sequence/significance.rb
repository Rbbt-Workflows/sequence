require 'rbbt/util/R'

module Sequence

  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input :exome, :boolean, "Limit analysis to exome bases", true
  dep do |jobname,options|
    if options[:exome]
      Sequence.job(:affected_genes, jobname, :organism => options[:organism], :mutations => options[:mutations])
    else
      Sequence.job(:genes, jobname, :organism => options[:organism], :positions => options[:mutations])
    end
  end
  task :binomial_significance => :tsv do |mutations, organism, exome|
    if Step === mutations
      Step.wait_for_jobs dependencies + [mutations]
      mutations = mutations.load 
    else
      Step.wait_for_jobs dependencies
    end

    if exome
      genes_step = step(:affected_genes)
      organism = genes_step.step(:mutated_isoforms_fast).inputs[:organism]
    else
      genes_step = step(:genes)
      organism = genes_step.inputs[:organism]
    end

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Matches", "Bases", "Frequency"], :namespace => organism, :type => :list, :cast => :to_f)

    mutation_count = Array === mutations ? 
      Misc.counts(mutations.collect{|m| m.split(":").values_at(0,1)*":"}) : 
      Misc.counts(Open.read(mutations).split("\n").collect{|m| m.split(":").values_at(0,1)*":"})

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

    mutations = mutations.to_a
    genes = genes.to_a

    num_mutations = mutation_count.values.inject(0){|acc,e| acc += e; acc}

    if exome
      total_bases = Organism.gene_list_exon_bases(genes)
      global_frequency = num_mutations.to_f / total_bases

      gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_exon_bases([gene]) }}
    else
      total_bases = Organism.gene_list_bases(genes)
      global_frequency = num_mutations.to_f / total_bases

      gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_bases([gene]) }}
    end

    TSV.traverse genes, :bar => "Computing significance" do |gene|
      mutations = gene_mutations[gene].uniq
      next if mutations.empty?
      matches = mutations.collect{|mutation| mutation_count[mutation]}.inject(0){|acc,e| acc += e}
      bases = gene2size[gene]
      next if bases == 0
      frequency = matches.to_f / bases

      tsv[gene] = [matches, bases, frequency]
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

    if exome
      total_bases = Organism.gene_list_bases(genes)
      global_frequency = num_mutations.to_f / total_bases

      gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_bases([gene]) }}
    else
      total_bases = Organism.gene_list_bases(genes)
      global_frequency = num_mutations.to_f / total_bases

      gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_bases([gene]) }}
    end

    TSV.traverse genes, :bar => "Computing significance" do |gene|
      mutations = gene_mutations[gene].uniq
      next if mutations.empty?
      matches = mutations.length
      bases = gene2size[gene]
      next if bases == 0
      frequency = matches.to_f / bases

      tsv[gene] = [matches, bases, frequency]
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

