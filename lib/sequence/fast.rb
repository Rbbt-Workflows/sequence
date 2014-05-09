module Sequence

  input *MUTATIONS_INPUT
  input *ORGANISM_INPUT
  input *WATSON_INPUT
  input *VCF_INPUT
  input *PRINCIPAL_INPUT
  dep &VCF_CONVERTER
  task :mutated_isoforms_fast => :tsv do |mutations,organism,watson,vcf,principal|
    if step(:genomic_mutations)
      Misc.consume_stream mutations, true
      mutations = step(:genomic_mutations)
    end

    transcript_protein = Sequence.transcript_protein(organism)
    exon_position = Sequence.exon_position(organism)
    exon_transcript_offsets = Sequence.exon_transcript_offsets(organism)

    chromosome_files = {}

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism
    dumper.init
    #log(:mutated_isoforms, :desc => "Mutated Isoforms") do |bar|

      TSV.traverse mutations, :cpus => 2, :bar => "Mutated Iso. Fast", :into => dumper, :type => :array do |mutation|
        next if mutation.nil?
        chr, pos, mut_str = mutation.split(":")
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
            offsets << [transcript, exon_offset.to_i + offset, strand] * ":"
          end if exon_transcript_offsets.include? exon
          offsets
        end.compact

        next if transcript_offsets.empty?

        mis = []
        mut_str.split(',').each do |mut|
          alleles = Sequence.alleles mut

          transcript_offsets.collect{|to| to.split ":" }.each do |transcript, transcript_offset, strand|
            next if principal and not Appris::PRINCIPAL_TRANSCRIPTS.include?(transcript)
            protein = transcript_protein[transcript]
            next if protein.nil? or protein.empty?

            begin
              codon = Sequence.codon_at_transcript_position(organism, transcript, transcript_offset.to_i);

              case codon

              when "UTR5", "UTR3"
                mis << [transcript, codon] * ":"

              else # Protein mutation
                triplet, offset, pos = codon.split ":"
                next if not triplet.length == 3
                original = Misc::CODON_TABLE[triplet]
                next if alleles.empty?
                pos = pos.to_i
                alleles.each do |allele|
                  change = case allele
                           when "Indel"
                             [original, pos + 1, "Indel"] * ""
                           when "FrameShift"
                             [original, pos + 1, "FrameShift"] * ""
                           else
                             allele = Misc::BASE2COMPLEMENT[allele] if watson and strand == "-1"
                             triplet[offset.to_i] = allele 
                             new = Misc::CODON_TABLE[triplet]
                             [original, pos + 1, new] * ""
                           end
                  mis << [protein, change] * ":"
                end
              end
            rescue TranscriptError
              Log.debug{$!.message}
            end
          end
        end
        next if mis.empty?

        [mutation, mis]
      end
    #end
  end
  export_synchronous :mutated_isoforms_fast

end
