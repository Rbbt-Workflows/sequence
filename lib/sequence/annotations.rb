module Sequence
  dep :mutated_isoforms_fast
  dep :splicing_mutations
  task :affected_genes => :tsv do
    organism = step(:mutated_isoforms_fast).inputs["organism"]

    protein_index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], :persist => true
    transcript_index = Organism.gene_transcripts(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism
    dumper.init
    TSV.traverse TSV.paste_streams([step(:splicing_mutations), step(:mutated_isoforms_fast)], :sort => true), :bar => "Affected Genes", :type => :array, :into => dumper do |line|
      next if line =~ /^#/
      m, *rest = line.split("\t")
      genes = rest.collect do |part|
        next if part.nil? or part.empty?
        e,c = part.split(":").first
        g = if e =~ /ENSP/
              next unless c =~ /^([A-Z*])\d+([A-Z*])/ and $1 != $2
              protein_index[e]
            else
              next if c =~ /^UTR\d$/
              transcript_index[e]
            end
        g
      end.compact.uniq

      next if genes.empty?
      [m, genes]
    end
  end
  export_asynchronous :affected_genes
end


