module Sequence

  dep do |jobname, options|
    options[:positions] = options[:mutations]
    Sequence.job(:reference_allele_at_genomic_positions, jobname, options)
  end
  dep do |jobname, options|
    options[:positions] = options[:mutations]
    Sequence.job(:genes_at_genomic_positions, jobname, options) if FalseClass === options[:watson]
  end
  input :mutations, :array, "Genomic Mutations"
  input :organism, :string, "Organism code", "Hsa"
  input :watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true
  task :mutation_reference => :tsv do |mutations,organism,watson|
    reference_allele_at_genomic_positions = step(:reference_allele_at_genomic_positions).load
    reference_allele_at_genomic_positions.unnamed = true

    if watson
      reference_allele_at_genomic_positions
    else
      new = {}
      gene_strand = Sequence.gene_strand_index(organism)
      mutation_genes = step(:genes_at_genomic_positions).load
      mutation_genes.unnamed = true
      log :fixing, "Fixing reference of crick strand muations"
      TSV.traverse reference_allele_at_genomic_positions, :into => new do |mutation,reference|
        base = mutation.split(":")[2]
        genes = mutation_genes[mutation]
        next if genes.nil? or genes.empty?

        reverse = genes.select{|gene| gene_strand[gene] == -1 }.any?
        forward = genes.select{|gene| gene_strand[gene] == 1 }.any?

        ref = case
              when (reverse and not forward)
                Misc::BASE2COMPLEMENT[reference]
              when (forward and not reverse)
                reference
              else
                base == reference ? Misc::BASE2COMPLEMENT[reference] : reference
              end
        [mutation,ref]
      end
      TSV.setup(new, :key_field => "Genomic Mutation", :fields => ["Mutation reference"], :type => :single)
    end
  end

  dep :mutation_reference
  input :mutations, :array, "Genomic Mutations"
  input :watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true
  task :mutation_type => :tsv do |mutations,watson|
    mutation_reference = step(:mutation_reference).load
    mutation_reference.unnamed = true

    mutation_type = {}
    TSV.traverse mutations, :into => mutation_type do |mutation|
      base = mutation.split(":")[2]
      reference = mutation_reference[mutation]

      type = case
             when (base.nil? or reference.nil? or base == "?" or reference == "?")
               "unknown"
             when base.index(',')
               "multiple"
             when base == reference
               "none"
             when (base.length > 1 or base == '-')
               "indel"
             when (not %w(A G T C).include? base and not %w(A G T C).include? reference) 
               "unknown"
             when (((Misc::IUPAC2BASE[base] || []) & ["A", "G"]).any? and     ((Misc::IUPAC2BASE[reference] || []) & ["T", "C"]).any?)
               "transversion"
             when (((Misc::IUPAC2BASE[base] || []) & ["T", "C"]).any? and     ((Misc::IUPAC2BASE[reference] || []) & ["A", "G"]).any?)
               "transversion"
             when (((Misc::IUPAC2BASE[base] || []) & ["A", "G"]).any? and     ((Misc::IUPAC2BASE[reference] || [nil]) & ["T", "C", nil]).empty?)
               "transition"
             when (((Misc::IUPAC2BASE[base] || []) & ["T", "C"]).any? and     ((Misc::IUPAC2BASE[reference] || [nil]) & ["A", "G", nil]).empty?)
               "transition"
             else
               "unknown"
             end
      [mutation, type]
    end

    TSV.setup(mutation_type, :key_field => "Genomic Mutation", :fields => ["Mutation type"], :type => :single)
  end
end
