module Sequence

  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input *WATSON_INPUT
  input *VCF_INPUT
  input *PRINCIPAL_INPUT
  input *NS_INPUT
  input *CODING_INPUT
  dep &VCF_CONVERTER
  task :mutated_isoforms_fast => :tsv do |mutations,organism,watson,vcf,principal,ns,coding|
    if dependencies.select{|d| d.task_name == :genomic_mutations}.any?
      mutations.close if IO === mutations
      mutations = step(:genomic_mutations)
    end

    raise ParameterException, "No mutations specified" if mutations.nil?

    blacklist_chromosomes = Organism.blacklist_chromosomes(organism).list
    transcript_protein = Sequence.transcript_protein(organism)
    exon_position = Sequence.exon_position(organism)
    exon_transcript_offsets = Sequence.exon_transcript_offsets(organism)
    biotype = Sequence.transcript_biotype(organism)

    chromosome_files = {}

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism
    dumper.init

    cpus = config('cpus', 'mutated_isoforms', :default => 3)
    TSV.traverse mutations, :cpus => cpus.to_i, :bar => self.progress_bar("Mutated Iso. Fast"), :into => dumper, :type => :array do |mutation|
      next if mutation.nil?
      chr, pos, mut_str = mutation.split(":")
      if blacklist_chromosomes.include? chr
        Log.low "Chromosome '#{chr}' is blacklisted"
        next
      end
      next if mut_str.nil?
      chr.sub!(/^chr/i,'')
      pos = pos.to_i

      index = chromosome_files[chr] ||= Sequence.exon_chromosome_index(organism, chr)
      exons = index[pos]

      next if exons.empty?

      transcript_offsets = exons.inject([]) do |offsets,exon|
        next offsets unless exon_position.include? exon
        strand, start, eend = exon_position[exon]
        if strand == 1
          offset = pos - start
        else
          offset = eend - pos
        end
        Misc.zip_fields(exon_transcript_offsets[exon]).each do |transcript, exon_offset|
          next if coding and biotype[transcript] != 'protein_coding'
          offsets << [transcript, exon_offset.to_i + offset, strand] * ":"
        end if exon_transcript_offsets.include? exon
        offsets
      end.compact

      next if transcript_offsets.empty?

      mis = []
      mut_str.split(/[,\/]/).each do |mut|
        alleles = Sequence.alleles mut

        transcript_offsets.collect{|to| to.split ":" }.each do |transcript, transcript_offset, strand|
          next if principal and organism =~ /^Hsa/ and not Appris::PRINCIPAL_TRANSCRIPTS.include?(transcript)
          protein = transcript_protein[transcript]
          next if protein.nil? or protein.empty?

          begin
            codon = Sequence.codon_at_transcript_position(organism, transcript, transcript_offset.to_i);

            case codon

            when "UTR5", "UTR3"
              mis << [transcript, codon] * ":"

            else # Protein mutation
              triplet, offset, pos = codon.split ":"
              if triplet.length < 3
                change = ["X",pos.to_i+1,"X"]*""
                mis << [protein, change] * ":"
              else
                original = Misc::CODON_TABLE[triplet] || 'X'
                next if alleles.empty?
                pos = pos.to_i
                alleles.each do |allele|
                  change = case allele
                          when "Indel"
                            change, pos, original = Sequence.downstream(organism, transcript, mut, codon)
                            [original, pos + 1, "Indel(#{change})"] * ""
                          when "FrameShift"
                            change, pos, original = Sequence.downstream(organism, transcript, mut, codon)
                            [original, pos + 1, "FrameShift(#{change})"] * ""
                          when /DNV\(([ATCG]+)\)/
                            as = $1.split("")
                            as = as.collect{|allele| Misc::BASE2COMPLEMENT[allele] } if watson and strand == "-1"
                            allele1, allele2 = as
                            offset1 = offset.to_i
                            offset2 = strand == "-1" ? offset1 - 1 : offset1 + 1
                            if offset2 < 0 or offset2 > 2 
                              change, pos, original = Sequence.downstream(organism, transcript, mut, codon)
                              [original, pos + 1, "Indel(#{change})"] * ""
                            else
                              triplet[offset1.to_i] = allele1
                              triplet[offset2.to_i] = allele2 
                              new = Misc::CODON_TABLE[triplet] || 'X'
                              [original, pos + 1, new] * ""
                            end
                          else
                            allele = Misc::BASE2COMPLEMENT[allele] if watson and strand == "-1"
                            triplet[offset.to_i] = allele 
                            new = Misc::CODON_TABLE[triplet] || 'X'
                            [original, pos + 1, new] * ""
                          end
                  mis << [protein, change] * ":"
                end
              end
            end
          rescue TranscriptError
            Log.debug{$!.message}
          end
        end
      end
      mis.reject!{|mi| mi !~ /ENS.*P\d+:([A-Z*]+)\d+([A-Z*]+)/i or $1 == $2 } if ns
      next if mis.empty?

      [mutation, mis]
    end
  end
  export_stream :mutated_isoforms_fast

end
