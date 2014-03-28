module Sequence

  dep :mutated_isoforms
  input :mutations, :array, "Genomic Mutations"
  input :organism, :string, "Organism code", "Hsa"
  task :non_synonymous => :tsv do |mutations,organism|
    tsv = {}
    gene2protein = Organism.transcripts(organism).index :target => "Ensembl Gene ID", :sources => "Ensembl Protein ID", :unnamed => true, :persist => true
    mutated_isoforms = step(:mutated_isoforms).load 
    mutated_isoforms.unnamed = true
    TSV.traverse mutations, :into => tsv do |mutation|
      mis = mutated_isoforms[mutation]
      next if mis.nil?
      values = []

      mutated_proteins = mis.collect{|mi| 
        protein, _sep, change = mi.partition ":"

        next unless change =~ /([A-Z*])\d+(.*)/ and $1 != $2

        protein
      }.compact

      next if mutated_proteins.empty?

      genes_with_mutated_isoforms = gene2protein.values_at(*mutated_proteins).compact.flatten.uniq
      values << genes_with_mutated_isoforms

      [mutation, values]
    end

    TSV.setup(tsv, :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :namespace => organism, :type => :flat) 
  end

  
  dep do |jobname, options|
    options[:positions] = options[:mutations]
    Sequence.job(:exon_junctions_at_genomic_positions, jobname, options)
  end
  dep :mutation_type
  input :mutations, :array, "Genomic Mutations"
  input :organism, :string, "Organism code", "Hsa"
  task :transcripts_with_affected_splicing => :tsv do |mutations, organism|
    return Transcript.setup([], "Ensembl Transcript ID", organism) if mutations.empty?

    exon_junctions_at_genomic_positions = step(:exon_junctions_at_genomic_positions).load
    exon_junctions_at_genomic_positions.unnamed = true
    mutation_type = step(:mutation_type).load
    mutation_type.unnamed = true

    exon_transcripts = Sequence.exon_transcripts(organism)
    transcript_exons = Sequence.transcript_exons(organism)

    tsv = {}
    TSV.traverse exon_junctions_at_genomic_positions, :into => tsv do |mutation,junctions|
      next if mutation_type[mutation] == "none" or junctions.empty?

      affected_transcripts = []
      junctions.each do |junction|
        exon, _sep, junction_type = junction.partition ":"
        transcripts = exon_transcripts[exon]
        next if transcripts.nil?

        transcripts.select! do |transcript|
          transcript_info = transcript_exons[transcript]
          next if transcript_info.nil?

          total_exons = transcript_info[0].length
          next if transcript_info[0].index(exon).nil?
          rank = transcript_info[1][transcript_info[0].index(exon)].to_i

          case
          when (rank == 1 and junction_type =~ /acceptor/)
            false
          when (rank == total_exons and junction_type =~ /donor/)
            false
          else
            true
          end
        end
        affected_transcripts.concat transcripts

      end
      next if  affected_transcripts.empty?
      [mutation, affected_transcripts]
    end
    TSV.setup(tsv, :key_field => "Genomic Mutation", :fields => ["Ensembl Transcript ID"], :namespace => organism, :type => :flat) 
  end

  dep :transcripts_with_affected_splicing
  dep :non_synonymous
  input :organism, :string, "Organism code", "Hsa"
  task :affected_genes => :tsv do |organism|
    transcripts_with_affected_splicing = step(:transcripts_with_affected_splicing).load
    non_synonymous = step(:non_synonymous).load
    transcripts_with_affected_splicing.unnamed = true
    non_synonymous.unnamed = true

    transcript_gene = Sequence.transcript_gene(organism)

    tsv = {}
    non_synonymous.each do |mutation,genes|
      next unless genes and genes.any?
      tsv[mutation] ||= []
      tsv[mutation].concat genes
    end

    transcripts_with_affected_splicing.each do |mutation,transcripts|
      next unless transcripts and transcripts.any?
      genes = transcript_gene.values_at(*transcripts).compact
      
      next if genes.empty?

      tsv[mutation] ||= []
      tsv[mutation].concat genes
    end

    tsv.each{|k,v| v.uniq! }

    TSV.setup(tsv, :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => organism)
  end

end
