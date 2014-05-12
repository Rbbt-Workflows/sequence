module Sequence
  dep :mutated_isoforms
  dep :splicing_mutations
  task :affected_genes => :tsv do
    organism = step(:mutated_isoforms).grace.inputs[1]

    protein_index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], :persist => true
    transcript_index = Organism.gene_transcripts(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse TSV.paste_streams([step(:mutated_isoforms), step(:splicing_mutations)], :sort => true), :bar => "Affected Genes", :type => :array, :into => dumper do |line|
      next if line =~ /^#/
      m, *rest = line.split("\t")
      genes = rest.collect do |part|
        e = part.split(":").first
        if e =~ /ENSP/
          protein_index[e]
        else
          transcript_index[e]
        end
      end.uniq
      next if genes.empty?
      [m, genes]
    end
  end
  export_asynchronous :affected_genes
end


