require 'rbbt/util/R'

module Sequence

  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input :exome, :boolean, "Limit analysis to exome bases", true
  input :num_samples, :integer, "Number of samples considered", 1
  dep do |jobname,options|
    if options[:exome]
      Sequence.job(:affected_genes, jobname, :organism => options[:organism], :mutations => options[:mutations])
    else
      Sequence.job(:genes, jobname, :organism => options[:organism], :positions => options[:mutations])
    end
  end
  task :binomial_significance => :tsv do |mutations, organism, exome, num_samples|
    if Step === mutations
      Step.wait_for_jobs dependencies + [mutations]
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

    genes = Set.new
    good_mutations = Set.new
    gene_positions = {}
    TSV.traverse genes_step, :type => :array, :bar => "Mutation catalog" do |mutation,g|
      next if mutation =~ /^#/
      mutation, *g = mutation.split("\t")
      next unless g and g.any?
      position = mutation.split(":").values_at(0,1)*":"
      good_mutations << position
      g.each do |gene|
        genes << gene
        gene_positions[gene] ||= []
        gene_positions[gene] << position
      end
    end

    position_counts = {}
    TSV.traverse mutations, :type => :array, :bar => "Selecting and counting matching mutations"  do |mutation|
      position = mutation.split(":").values_at(0,1)*":"
      next unless good_mutations.include? position
      position_counts[position] ||= 0
      position_counts[position] += 1
    end

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Matches", "Bases", "Frequency"], :namespace => organism, :type => :list, :cast => :to_f)

    genes = genes.to_a

    num_mutations = position_counts.values.inject(0){|acc,e| acc += e; acc}

    log :computing_sizes
    if exome
      total_bases = Organism.gene_list_exon_bases(genes) * num_samples
      global_frequency = num_mutations.to_f / total_bases

      gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_exon_bases([gene]) * num_samples }}
    else
      total_bases = Organism.gene_list_bases(genes) * num_samples
      global_frequency = num_mutations.to_f / total_bases

      gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_bases([gene]) * num_samples }}
    end

    TSV.traverse genes, :bar => "Computing significance", :type => :array do |gene|
      positions = gene_positions[gene].uniq
      next if positions.empty?
      matches = positions.collect{|position| position_counts[position]}.compact.inject(0){|acc,e| acc += e}
      next if matches < 2

      bases = gene2size[gene]
      next if bases == 0
      frequency = matches.to_f / bases
      if matches >= bases
        log :skip, "Gene #{ gene } has more matches than mutations. Skipping. Are you sure these mutations come from #{num_samples} samples?"
        next
      end

      tsv[gene] = [matches, bases, frequency]
    end

    begin
      new = tsv.R "
      data = cbind(data,p.value=apply(data, 1, function(v){v = as.numeric(v); binom.test(v[1], v[2], #{ global_frequency }, 'greater')$p.value}))
      ", :key => "Ensembl Gene ID" 
      new.namespace = organism
      new
    rescue Exception
      Log.exception $!
      tsv.annotate({}).add_field "p.value" do
      end
    end
  end

  #input *MUTATIONS_INPUT
  #input *ORGANISM_INPUT
  #input :exome, :boolean, "Limit analysis to exome bases", true
  #dep do |jobname,options|
  #  if options[:exome]
  #    Sequence.job(:affected_genes, jobname, :organism => options[:organism], :mutations => options[:mutations])
  #  else
  #    Sequence.job(:genes, jobname, :organism => options[:organism], :positions => options[:mutations])
  #  end
  #end
  #task :binomial_significance_old => :tsv do |mutations, organism, exome|
  #  if Step === mutations
  #    Step.wait_for_jobs dependencies + [mutations]
  #    log :loading_mutations
  #    mutations = mutations.load 
  #  else
  #    Step.wait_for_jobs dependencies
  #  end

  #  log :get_feature_job
  #  if exome
  #    genes_step = step(:affected_genes)
  #    organism = genes_step.step(:mutated_isoforms_fast).inputs[:organism]
  #  else
  #    genes_step = step(:genes)
  #    organism = genes_step.inputs[:organism]
  #  end

  #  log :getting_matched_mutations
  #  matched_mutations = TSV.traverse genes_step, :type => :array, :into => [], :bar => true do |line|
  #    mutation, genes = line.partition("\t")
  #    next unless genes
  #    mutation
  #  end

  #  matched_mutations.uniq!

  #  tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Matches", "Bases", "Frequency"], :namespace => organism, :type => :list, :cast => :to_f)

  #  good_mutations = Array === mutations ? 
  #    (mutations & matched_mutations) :
  #    (Open.read(mutations).split("\n") & matched_mutations)

  #  positions = good_mutations.collect{|m| m.split(":").values_at(0,1)*":"} 

  #  mutation_count = Misc.counts(positions)

  #  log :mutation_catalog
  #  genes = Set.new
  #  mutations = Set.new
  #  gene_mutations = {}
  #  TSV.traverse genes_step do |position,g|
  #    position = position.split(":").values_at(0,1)*":"
  #    mutations << position
  #    next if g.nil?
  #    g.each do |gene|
  #      genes << gene
  #      gene_mutations[gene] ||= []
  #      gene_mutations[gene] << position
  #    end
  #  end

  #  mutations = mutations.to_a
  #  genes = genes.to_a

  #  num_mutations = mutation_count.values.inject(0){|acc,e| acc += e; acc}

  #  if exome
  #    total_bases = Organism.gene_list_exon_bases(genes)
  #    global_frequency = num_mutations.to_f / total_bases

  #    gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_exon_bases([gene]) }}
  #  else
  #    total_bases = Organism.gene_list_bases(genes)
  #    global_frequency = num_mutations.to_f / total_bases

  #    gene2size = Misc.process_to_hash(genes){|genes| genes.collect{|gene| Organism.gene_list_bases([gene]) }}
  #  end

  #  TSV.traverse genes, :bar => "Computing significance" do |gene|
  #    mutations = gene_mutations[gene].uniq
  #    next if mutations.empty?
  #    matches = mutations.collect{|mutation| mutation_count[mutation]}.inject(0){|acc,e| acc += e}
  #    bases = gene2size[gene]
  #    next if bases == 0
  #    frequency = matches.to_f / bases

  #    tsv[gene] = [matches, bases, frequency]
  #  end

  #  begin
  #    new = tsv.R "
  #    data = cbind(data,p.value=apply(data, 1, function(v){v = as.numeric(v); binom.test(v[1], v[2], #{ global_frequency }, 'greater')$p.value}))
  #    ", :key => "Ensembl Gene ID" 
  #    new.namespace = organism
  #    new
  #  rescue Exception
  #    tsv.annotate({}).add_field "p.value" do
  #    end
  #  end
  #end
  #export_asynchronous :binomial_significance

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

