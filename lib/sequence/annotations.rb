module Sequence
  dep :mutated_isoforms
  task :affected_genes => :tsv do
    organism = step(:mutated_isoforms).inputs[1]

    index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], :persist => true

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse step(:mutated_isoforms), :into => dumper do |m,mis|
      next if mis.nil? or mis.empty?
      proteins = mis.collect{|m| m.partition(":")[0]}.select{|p| p =~ /^ENSP/ }.uniq
      genes = index.values_at(*proteins).compact
      next if genes.empty?
      [m, genes.uniq]
    end
  end
  export_asynchronous :affected_genes
end


